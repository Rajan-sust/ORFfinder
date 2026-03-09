package main

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"
)

// ─────────────────────────────────────────────
//  Data Structures
// ─────────────────────────────────────────────

// ORF holds a single open reading frame result.
type ORF struct {
	SeqID    string
	Strand   string // "+" or "-"
	Frame    int    // 0, 1, or 2
	Start    int    // 1-based, inclusive
	End      int    // 1-based, inclusive
	Length   int    // nucleotide length
	Protein  string // translated amino-acid sequence (optional)
	NucSeq   string // raw nucleotide sequence (optional)
}

// searchJob bundles a single sequence for worker goroutines.
type searchJob struct {
	id  string
	seq []byte
}

// ─────────────────────────────────────────────────────────────────
//  Genetic Code — NCBI Translation Table 11
//  Bacteria, Archaea, Prokaryotic Viruses & Chloroplast Plastids
//  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
//
//  Encodes the ELONGATION meaning of every codon (64 total).
//  Start-codon recognition (ATG, GTG, TTG, ATT, ATC, ATA)
//  is handled separately in isStart() / scanStrand().
// ─────────────────────────────────────────────────────────────────

var codonTable = map[string]byte{
	// ── T-second-base column ────────────────────────────────────
	"TTT": 'F', "TTC": 'F', // Phe
	"TTA": 'L', "TTG": 'L', // Leu
	"TCT": 'S', "TCC": 'S', "TCA": 'S', "TCG": 'S', // Ser
	"TAT": 'Y', "TAC": 'Y', // Tyr
	"TAA": '*', "TAG": '*', // Ter — stop
	"TGT": 'C', "TGC": 'C', // Cys
	"TGA": '*',               // Ter — stop
	"TGG": 'W',               // Trp

	// ── C-second-base column ────────────────────────────────────
	"CTT": 'L', "CTC": 'L', "CTA": 'L', "CTG": 'L', // Leu
	"CCT": 'P', "CCC": 'P', "CCA": 'P', "CCG": 'P', // Pro
	"CAT": 'H', "CAC": 'H', // His
	"CAA": 'Q', "CAG": 'Q', // Gln
	"CGT": 'R', "CGC": 'R', "CGA": 'R', "CGG": 'R', // Arg

	// ── A-second-base column ────────────────────────────────────
	"ATT": 'I', "ATC": 'I', "ATA": 'I', // Ile  (also alt start in Table 11)
	"ATG": 'M',               // Met  — primary start codon
	"ACT": 'T', "ACC": 'T', "ACA": 'T', "ACG": 'T', // Thr
	"AAT": 'N', "AAC": 'N', // Asn
	"AAA": 'K', "AAG": 'K', // Lys
	"AGT": 'S', "AGC": 'S', // Ser
	"AGA": 'R', "AGG": 'R', // Arg

	// ── G-second-base column ────────────────────────────────────
	"GTT": 'V', "GTC": 'V', "GTA": 'V', "GTG": 'V', // Val  (also alt start in Table 11)
	"GCT": 'A', "GCC": 'A', "GCA": 'A', "GCG": 'A', // Ala
	"GAT": 'D', "GAC": 'D', // Asp
	"GAA": 'E', "GAG": 'E', // Glu
	"GGT": 'G', "GGC": 'G', "GGA": 'G', "GGG": 'G', // Gly
}

// isStop returns true for the three universal stop codons.
func isStop(c []byte) bool {
	s := string(c)
	return s == "TAA" || s == "TAG" || s == "TGA"
}

// isStart returns true for ATG (canonical) only.
// Table 11 allows GTG / TTG / ATT / ATC / ATA as alternative starts,
// but only when they are the FIRST codon of a CDS. Since we scan
// for all ATG-initiated ORFs per standard ORFfinder behaviour, we
// keep ATG-only here. Pass -altstart flag (future) to enable the rest.
func isStart(c []byte) bool {
	s := string(c)
	return s == "ATG"
}

// translate converts a nucleotide slice to an amino-acid string.
// Assumes seq starts at codon boundary; stops at first stop codon.
func translate(seq []byte) string {
	var sb strings.Builder
	sb.Grow(len(seq) / 3)
	for i := 0; i+2 < len(seq); i += 3 {
		codon := string(seq[i : i+3])
		aa, ok := codonTable[codon]
		if !ok {
			aa = 'X'
		}
		if aa == '*' {
			break
		}
		sb.WriteByte(aa)
	}
	return sb.String()
}

// ─────────────────────────────────────────────
//  Reverse Complement
// ─────────────────────────────────────────────

var complement = [256]byte{
	'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
	'a': 'T', 't': 'A', 'g': 'C', 'c': 'G',
	'N': 'N', 'n': 'N',
}

// revComp returns the reverse complement of a DNA sequence.
func revComp(seq []byte) []byte {
	n := len(seq)
	rc := make([]byte, n)
	for i, b := range seq {
		c := complement[b]
		if c == 0 {
			c = 'N'
		}
		rc[n-1-i] = c
	}
	return rc
}

// ─────────────────────────────────────────────
//  Core ORF Scanner
// ─────────────────────────────────────────────

// scanStrand finds all ORFs in a single strand (already in 5'→3' orientation).
// strand is "+" or "-"; offset is used to re-map "-" strand positions back to
// forward-strand coordinates.
// For each stop codon only the LONGEST ORF (earliest ATG) is reported.
func scanStrand(
	id string,
	seq []byte,
	strand string,
	seqLen int,
	minLen int,
	translate bool,
) []ORF {
	results := make([]ORF, 0, 64)

	for frame := 0; frame < 3; frame++ {
		// firstStart tracks the earliest ATG seen since the last stop codon.
		// -1 means no ATG has been seen yet in this open window.
		firstStart := -1

		for pos := frame; pos+2 < len(seq); pos += 3 {
			codon := seq[pos : pos+3]
			if isStart(codon) {
				if firstStart == -1 {
					firstStart = pos // record only the FIRST (longest) ATG
				}
			} else if isStop(codon) {
				if firstStart != -1 {
					stopEnd := pos + 3 // exclusive end of stop codon
					nucLen := stopEnd - firstStart
					if nucLen >= minLen {
						orf := ORF{
							SeqID:  id,
							Strand: strand,
							Frame:  frame + 1,
							Length: nucLen,
						}
						nucSeq := seq[firstStart:stopEnd]

						if strand == "+" {
							orf.Start = firstStart + 1
							orf.End = stopEnd
						} else {
							// Map back to forward-strand 1-based coordinates
							orf.Start = seqLen - stopEnd + 1
							orf.End = seqLen - firstStart
						}

						if translate {
							orf.Protein = translateSeq(nucSeq)
							orf.NucSeq = string(nucSeq)
						}
						results = append(results, orf)
					}
				}
				// Reset window — any ATG after this stop starts a new ORF.
				firstStart = -1
			}
		}
	}
	return results
}

// translateSeq translates including the stop codon (excluded from protein).
func translateSeq(seq []byte) string {
	var sb strings.Builder
	sb.Grow(len(seq) / 3)
	for i := 0; i+2 < len(seq); i += 3 {
		codon := string(seq[i : i+3])
		aa, ok := codonTable[codon]
		if !ok {
			aa = 'X'
		}
		if aa == '*' {
			break
		}
		sb.WriteByte(aa)
	}
	return sb.String()
}

// findORFs searches both strands of a sequence concurrently.
func findORFs(id string, seq []byte, minLen int, doTranslate bool) []ORF {
	seqLen := len(seq)

	var wg sync.WaitGroup
	ch := make(chan []ORF, 2)

	// Forward strand
	wg.Add(1)
	go func() {
		defer wg.Done()
		ch <- scanStrand(id, seq, "+", seqLen, minLen, doTranslate)
	}()

	// Reverse strand
	wg.Add(1)
	go func() {
		defer wg.Done()
		rc := revComp(seq)
		ch <- scanStrand(id, rc, "-", seqLen, minLen, doTranslate)
	}()

	go func() {
		wg.Wait()
		close(ch)
	}()

	var all []ORF
	for orfs := range ch {
		all = append(all, orfs...)
	}
	return all
}

// ─────────────────────────────────────────────
//  FASTA Parser (streaming, zero-copy-friendly)
// ─────────────────────────────────────────────

type FastaRecord struct {
	ID  string
	Seq []byte
}

// parseFasta reads FASTA records from a file path and sends them to the
// returned channel. Parsing is done in a dedicated goroutine.
func parseFasta(path string) (<-chan FastaRecord, <-chan error) {
	out := make(chan FastaRecord, 32)
	errc := make(chan error, 1)

	go func() {
		defer close(out)
		defer close(errc)

		f, err := os.Open(path)
		if err != nil {
			errc <- err
			return
		}
		defer f.Close()

		var currentID string
		var seqBuf []byte

		flush := func() {
			if currentID != "" && len(seqBuf) > 0 {
				// Uppercase and validate
				upper := make([]byte, len(seqBuf))
				for i, b := range seqBuf {
					if b >= 'a' && b <= 'z' {
						upper[i] = b - 32
					} else {
						upper[i] = b
					}
				}
				out <- FastaRecord{ID: currentID, Seq: upper}
			}
		}

		scanner := bufio.NewScanner(f)
		scanner.Buffer(make([]byte, 8*1024*1024), 8*1024*1024) // 8 MB line buffer

		for scanner.Scan() {
			line := scanner.Text()
			if len(line) == 0 {
				continue
			}
			if line[0] == '>' {
				flush()
				seqBuf = seqBuf[:0] // reset, reuse backing array if possible
				// Take only the first token as ID
				fields := strings.Fields(line[1:])
				if len(fields) > 0 {
					currentID = fields[0]
				} else {
					currentID = "unknown"
				}
			} else {
				seqBuf = append(seqBuf, []byte(line)...)
			}
		}
		flush()

		if err := scanner.Err(); err != nil {
			errc <- err
		}
	}()

	return out, errc
}

// ─────────────────────────────────────────────
//  Output Writers
// ─────────────────────────────────────────────

// writeGFF3 writes ORF results in GFF3 format.
func writeGFF3(w *bufio.Writer, orfs []ORF) {
	fmt.Fprintln(w, "##gff-version 3")
	for i, o := range orfs {
		fmt.Fprintf(w, "%s\tORFfinder\tCDS\t%d\t%d\t.\t%s\t%d\tID=orf%d;length=%d\n",
			o.SeqID, o.Start, o.End, o.Strand, o.Frame-1, i+1, o.Length)
	}
}

// writeTSV writes ORF results as tab-separated values.
func writeTSV(w *bufio.Writer, orfs []ORF, showSeqs bool) {
	if showSeqs {
		fmt.Fprintln(w, "SeqID\tStrand\tFrame\tStart\tEnd\tLength\tProtein\tNucleotide")
	} else {
		fmt.Fprintln(w, "SeqID\tStrand\tFrame\tStart\tEnd\tLength")
	}
	for _, o := range orfs {
		if showSeqs {
			fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n",
				o.SeqID, o.Strand, o.Frame, o.Start, o.End, o.Length, o.Protein, o.NucSeq)
		} else {
			fmt.Fprintf(w, "%s\t%s\t%d\t%d\t%d\t%d\n",
				o.SeqID, o.Strand, o.Frame, o.Start, o.End, o.Length)
		}
	}
}

// writeFasta writes translated protein sequences in FASTA format.
func writeFasta(w *bufio.Writer, orfs []ORF) {
	for i, o := range orfs {
		fmt.Fprintf(w, ">orf%d seqid=%s strand=%s frame=%d start=%d end=%d len=%d\n",
			i+1, o.SeqID, o.Strand, o.Frame, o.Start, o.End, o.Length)
		// Wrap at 60 characters
		p := o.Protein
		for len(p) > 60 {
			fmt.Fprintln(w, p[:60])
			p = p[60:]
		}
		if len(p) > 0 {
			fmt.Fprintln(w, p)
		}
	}
}

// ─────────────────────────────────────────────
//  Overlap Filter
// ─────────────────────────────────────────────

// removeOverlaps filters ORFs using strand-aware greedy interval scheduling.
//
// Three-tier biological rules:
//
//  TIER 1 — OPPOSITE-STRAND overlaps → ALWAYS allowed.
//    Antisense overlapping genes are a well-documented prokaryotic feature.
//    Two genes on opposite strands are transcribed and translated independently
//    even when they share a genomic coordinate range.
//
//  TIER 2 — SAME-STRAND CROSSING overlaps within -maxoverlap → ALLOWED.
//    "Crossing" means the intervals interleave at their boundary:
//         gene A:  ──────────[endA]
//         gene B:        [startB]──────────
//    This is translational coupling: the start of gene B overlaps the end
//    of gene A by a few nucleotides (classically 1–4 nt, e.g. the "ATGA"
//    motif where TGA stop of gene A and ATG start of gene B share one base).
//    The -maxoverlap flag (default 4) sets the maximum crossing nt allowed.
//    Pseudomonas aeruginosa example:
//         gene 1  3169..4278
//         gene 2  4275..6695  ← 4-nt crossing overlap → translational coupling ✓
//
//  TIER 3 — all other same-strand overlaps → CONFLICT, keep longest.
//    This includes:
//      • Nested ORFs (one fully or largely inside the other) — not coupling.
//      • Crossing overlaps larger than -maxoverlap — not coupling.
//    The longer ORF wins; the shorter is discarded.
//
// Crossing vs nesting — the key distinction:
//    Crossing:  a.Start < b.Start < a.End < b.End  (or mirror)  → may be coupling
//    Nesting:   a.Start < b.Start AND b.End < a.End (or mirror)  → always conflict
//
// Algorithm — greedy interval scheduling by length (longest-first):
//   1. Sort ORFs by length descending; tie-break by Start ascending.
//   2. For each candidate ORF, check every same-strand accepted interval:
//        a. Compute overlap length = max(0, min(endA,endB) - max(startA,startB) + 1)
//        b. If overlap == 0 → no conflict.
//        c. If intervals CROSS and overlap <= maxSameStrandOverlap → coupling, no conflict.
//        d. Otherwise → conflict, discard candidate.
//   3. Re-sort accepted set by SeqID + Start for output.
func removeOverlaps(orfs []ORF, maxSameStrandOverlap int) []ORF {
	// Sort by length descending; tie-break by Start ascending (deterministic).
	sort.Slice(orfs, func(i, j int) bool {
		if orfs[i].Length != orfs[j].Length {
			return orfs[i].Length > orfs[j].Length
		}
		return orfs[i].Start < orfs[j].Start
	})

	type interval struct{ start, end int }
	accepted := make(map[string][]interval, 32) // key = seqid + "\x00" + strand
	kept := make([]ORF, 0, len(orfs))

	for _, o := range orfs {
		// Tier 1: opposite-strand ORFs go into a separate bucket → never compete.
		key := o.SeqID + "\x00" + o.Strand
		conflict := false
		for _, iv := range accepted[key] {
			overlap := min(o.End, iv.end) - max(o.Start, iv.start) + 1
			if overlap <= 0 {
				continue // no overlap at all — fine
			}
			// Tier 2: crossing overlap within coupling tolerance → fine.
			// Crossing means the intervals interleave at their ends, NOT nest.
			crossing := (o.Start < iv.start && iv.start <= o.End && o.End < iv.end) ||
				(iv.start < o.Start && o.Start <= iv.end && iv.end < o.End)
			if crossing && overlap <= maxSameStrandOverlap {
				continue // translational coupling — allow both
			}
			// Tier 3: any other same-strand overlap is a genuine conflict.
			conflict = true
			break
		}
		if !conflict {
			accepted[key] = append(accepted[key], interval{o.Start, o.End})
			kept = append(kept, o)
		}
	}

	// Re-sort by SeqID + Start for final output.
	sort.Slice(kept, func(i, j int) bool {
		if kept[i].SeqID != kept[j].SeqID {
			return kept[i].SeqID < kept[j].SeqID
		}
		return kept[i].Start < kept[j].Start
	})
	return kept
}

// min/max helpers (Go 1.21+ has these built-in; keep for compatibility).
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// ─────────────────────────────────────────────
//  Pipeline Orchestrator
// ─────────────────────────────────────────────

func run(
	inputPath string,
	outputPath string,
	minLen int,
	outFormat string,
	doTranslate bool,
	numWorkers int,
	maxOverlap int,
	verbose bool,
) error {
	start := time.Now()

	// ── Parse FASTA (producer) ──
	records, parseErrs := parseFasta(inputPath)

	// ── Worker pool (consumers) ──
	jobs := make(chan searchJob, numWorkers*2)
	resultCh := make(chan []ORF, numWorkers*2)
	var workerWg sync.WaitGroup

	for i := 0; i < numWorkers; i++ {
		workerWg.Add(1)
		go func() {
			defer workerWg.Done()
			for job := range jobs {
				orfs := findORFs(job.id, job.seq, minLen, doTranslate)
				if len(orfs) > 0 {
					resultCh <- orfs
				}
			}
		}()
	}

	// ── Feed jobs from FASTA channel ──
	var feedWg sync.WaitGroup
	feedWg.Add(1)
	go func() {
		defer feedWg.Done()
		defer close(jobs)
		seqCount := 0
		for rec := range records {
			seqCount++
			jobs <- searchJob{id: rec.ID, seq: rec.Seq}
		}
		if verbose {
			fmt.Fprintf(os.Stderr, "[info] Dispatched %d sequences to %d workers\n", seqCount, numWorkers)
		}
	}()

	// ── Close resultCh when all workers finish ──
	go func() {
		feedWg.Wait()
		workerWg.Wait()
		close(resultCh)
	}()

	// ── Collect results ──
	var allORFs []ORF
	for batch := range resultCh {
		allORFs = append(allORFs, batch...)
	}

	// Check parse errors
	if err := <-parseErrs; err != nil {
		return fmt.Errorf("FASTA parse error: %w", err)
	}

	// ── Sort by SeqID, then Start position ──
	sort.Slice(allORFs, func(i, j int) bool {
		if allORFs[i].SeqID != allORFs[j].SeqID {
			return allORFs[i].SeqID < allORFs[j].SeqID
		}
		return allORFs[i].Start < allORFs[j].Start
	})

	// ── Remove overlapping ORFs (keep longest per locus) ──
	beforeFilter := len(allORFs)
	allORFs = removeOverlaps(allORFs, maxOverlap)

	elapsed := time.Since(start)

	if verbose {
		fmt.Fprintf(os.Stderr, "[info] Found %d ORFs before overlap filter, %d after (%d removed) in %s\n",
			beforeFilter, len(allORFs), beforeFilter-len(allORFs), elapsed)
	}

	// ── Write output ──
	var outFile *os.File
	var err error
	if outputPath == "-" || outputPath == "" {
		outFile = os.Stdout
	} else {
		outFile, err = os.Create(outputPath)
		if err != nil {
			return fmt.Errorf("cannot create output file: %w", err)
		}
		defer outFile.Close()
	}

	bw := bufio.NewWriterSize(outFile, 1<<20) // 1 MB write buffer
	defer bw.Flush()

	switch strings.ToLower(outFormat) {
	case "gff", "gff3":
		writeGFF3(bw, allORFs)
	case "fasta", "fa", "prot":
		writeFasta(bw, allORFs)
	case "tsv":
		writeTSV(bw, allORFs, doTranslate)
	default:
		writeTSV(bw, allORFs, doTranslate)
	}

	if verbose && outputPath != "-" && outputPath != "" {
		fmt.Fprintf(os.Stderr, "[info] Results written to %s\n", outputPath)
	}

	return nil
}

// ─────────────────────────────────────────────
//  Main / CLI
// ─────────────────────────────────────────────

func main() {
	var (
		inputPath   = flag.String("i", "", "Input FASTA file (required)")
		outputPath  = flag.String("o", "-", "Output file (default: stdout)")
		minLen      = flag.Int("min", 300, "Minimum ORF length in nucleotides (default: 300)")
		outFormat   = flag.String("fmt", "tsv", "Output format: tsv | gff3 | fasta  (default: tsv)")
		doTranslate = flag.Bool("translate", false, "Translate ORFs to protein sequences (slower)")
		numWorkers      = flag.Int("workers", runtime.NumCPU(), "Number of parallel worker goroutines")
		maxOverlap      = flag.Int("maxoverlap", 4, "Max same-strand overlap (nt) allowed between adjacent genes (default: 4, covers translational coupling)")
		verbose         = flag.Bool("v", false, "Verbose / progress output to stderr")
	)
	flag.Usage = func() {
		fmt.Fprintln(os.Stderr, `
╔══════════════════════════════════════════════════╗
║         ORFfinder — Prokaryotic ORF Scanner       ║
║         Concurrent Go implementation              ║
╚══════════════════════════════════════════════════╝

Usage:
  orffinder -i <input.fasta> [options]

Options:`)
		flag.PrintDefaults()
		fmt.Fprintln(os.Stderr, `
Examples:
  orffinder -i genome.fna -min 300 -fmt tsv -v
  orffinder -i genome.fna -min 150 -fmt gff3 -o orfs.gff
  orffinder -i genome.fna -translate -fmt fasta -o proteins.fa
  orffinder -i genome.fna -workers 16 -v
`)
	}
	flag.Parse()

	if *inputPath == "" {
		fmt.Fprintln(os.Stderr, "[error] Input file required (-i flag)")
		flag.Usage()
		os.Exit(1)
	}
	if *minLen < 3 {
		fmt.Fprintln(os.Stderr, "[error] Minimum ORF length must be >= 3 nt")
		os.Exit(1)
	}
	if *numWorkers < 1 {
		*numWorkers = 1
	}

	if *verbose {
		fmt.Fprintf(os.Stderr, "[info] ORFfinder starting — workers=%d minLen=%d format=%s maxoverlap=%d\n",
			*numWorkers, *minLen, *outFormat, *maxOverlap)
	}

	if err := run(*inputPath, *outputPath, *minLen, *outFormat, *doTranslate, *numWorkers, *maxOverlap, *verbose); err != nil {
		fmt.Fprintf(os.Stderr, "[error] %v\n", err)
		os.Exit(1)
	}
}