package main

import (
	"fmt"

	"github.com/biogo/biogo/align"
	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/seq/linear"
)

// struct for input reads

// Output struct for processed reads

// Slice struct to hold the slice
type Slice struct {
	X int
	Y int
}

func ProcessRead(readInfo ReadHolder) (Trios, error) {
	// get the read sequence
	read1 := readInfo.Seq1
	read2 := readInfo.Seq2
	// get the read quality

	// Initialize the output struct
	var output Trios
	// First get the cell barcodes and
	// Call fuzzy match function to locate call the barcodes
	// get the cell barcode
	cellBC := read1[0:16]
	// get the UMI
	umi := read1[16:28]
	// Check the Q30 score of the read
	output.CellBC = cellBC
	output.Umi = umi
	// Fuzzy match for thr tripBC
	beforePBC := "AAGTAATCTAGA"
	afterPBC := "GTCGAGATAA"
	beforeRBC := "CTATACGAAGTTATG"
	afterRBC := "GCTTTAAGGCCGGTCC"
	coordinates, err := FuzzyMatch(beforePBC, afterPBC, beforeRBC, afterRBC, read2)
	if err != nil {
		output := Trios{}
		return output, err
	}
	// et the TBC
	tbc := read2[coordinates[0]:coordinates[1]]
	rbc := read2[coordinates[2]:coordinates[3]]
	output.RBC = rbc
	output.TBC = tbc
	// create the output struct
	return output, nil
}

func FuzzyMatch(beforePBC, afterPBC, beforeRBC, afterRBC, read string) ([4]int8, error) {
	// Use sequence to find the pBC in the reads
	fsa := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(read))}
	fsa.Alpha = alphabet.DNAgapped
	fsb := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(beforePBC))}
	fsb.Alpha = alphabet.DNAgapped
	fsc := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(afterPBC))}
	fsc.Alpha = alphabet.DNAgapped
	// Use sequence to find the rBC in the reads
	fsd := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(beforeRBC))}
	fsd.Alpha = alphabet.DNAgapped
	fse := &linear.Seq{Seq: alphabet.BytesToLetters([]byte(afterRBC))}
	fse.Alpha = alphabet.DNAgapped
	// Here we penalize gaps (does not happen in sequencing errors) and we
	// like to penalize mismatches.
	fitted := align.Fitted{
		{0, -100, -100, -100, -100},
		{-100, 100, -10, -10, -10},
		{-100, -10, 100, -10, -10},
		{-100, -10, -10, 100, -10},
		{-100, -10, -10, -10, 100},
	}
	aln_1, err := fitted.Align(fsa, fsb)
	// Check if the alignment worked, if it returns error, print error
	if err != nil {
		empty := [4]int8{0, 0, 0, 0}
		return empty, err
	}
	// If the length of alignment is more than 1, discard the read
	if len(aln_1) > 1 {
		return [4]int8{0, 0, 0, 0}, fmt.Errorf("alignment length is more than 1")
	}
	aln_2, err := fitted.Align(fsa, fsc)
	if err != nil {
		empty := [4]int8{0, 0, 0, 0}
		return empty, err
	}
	// If the length of alignment is more than 1, discard the read
	if len(aln_2) > 1 {
		return [4]int8{0, 0, 0, 0}, fmt.Errorf("alignment length is more than 1")
	}
	// align the reads around rBC
	aln_3, err := fitted.Align(fsa, fsd)
	if err != nil {
		empty := [4]int8{0, 0, 0, 0}
		return empty, err
	}
	if len(aln_3) > 1 {
		return [4]int8{0, 0, 0, 0}, fmt.Errorf("alignment length is more than 1 for rBC")
	}
	aln_4, err := fitted.Align(fsa, fse)
	if err != nil {
		empty := [4]int8{0, 0, 0, 0}
		return empty, err
	}
	if len(aln_4) > 1 {
		return [4]int8{0, 0, 0, 0}, fmt.Errorf("alignment length is more than 1 for rBC")
	}
	// If the length of alignment is 1, return the coordinates
	var coordinates [4]int8
	coordinates[0] = int8(aln_1[0].Features()[0].End())
	coordinates[1] = int8(aln_2[0].Features()[0].Start())
	coordinates[2] = int8(aln_3[0].Features()[0].End())
	coordinates[3] = int8(aln_4[0].Features()[0].Start())
	// Check if the length of those two are correct
	length := coordinates[1] - coordinates[0]
	if length != 12 {
		return [4]int8{0, 0, 0, 0}, fmt.Errorf("wrong length for pBC")
	}
	// Check if the length of the rBC is correct
	length_rBC := coordinates[3] - coordinates[2]
	if length_rBC != 25 {
		return [4]int8{0, 0, 0, 0}, fmt.Errorf("wrong length for rBC")
	}
	// Return the coordinates of the aligned sequence
	return coordinates, nil
}
