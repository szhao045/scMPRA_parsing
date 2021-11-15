package main

import (
	"bufio"
	"compress/gzip"
	"fmt"
	"os"
	"runtime"
)

type ReadHolder struct {
	Seq1 string
	Seq2 string
}
type PairedRead struct {
	R1 string
	R2 string
}

type Trios struct {
	CellBC string
	Umi    string
	TBC    string
	RBC    string
}

func reader(pair PairedRead) map[Trios]int {
	read1 := pair.R1
	read2 := pair.R2
	fmt.Println("Reading Read 1 file:", read1)
	fmt.Println("Reading Read 2 file:", read2)
	// Open the file
	f1, err := os.Open(read1)
	if err != nil {
		fmt.Println(err)
		return nil
	}
	// Try to unzip the file
	f1_unzip, err := gzip.NewReader(f1)
	if err != nil {
		fmt.Println("The error is", err)
		return nil
	}
	f2, err := os.Open(read2)
	if err != nil {
		fmt.Println(err)
		return nil
	}
	f2_unzip, err := gzip.NewReader(f2)
	if err != nil {
		fmt.Println(err)
		return nil
	}
	// Close the file and garbage collect the file handle
	defer f1.Close()
	defer f2.Close()
	scanner1 := bufio.NewScanner(f1_unzip)
	scanner2 := bufio.NewScanner(f2_unzip)
	// Initialize a counter
	line := 0
	wrong_trio_counter := 0
	// Initialize a map to hold the data
	trio_holder := make(map[Trios]int)
	// Scan through the file line by line
	for scanner1.Scan() && scanner2.Scan() {
		line++
		// Initiate a holder of read sequences and qualities
		var readInfo ReadHolder
		// Add
		if line%4 == 2 {
			read1 := scanner1.Text()
			readInfo.Seq1 = read1
			read2 := scanner2.Text()
			readInfo.Seq2 = read2

			// Pass readinfo to the function to process the read
			trio, err := ProcessRead(readInfo)
			// Handles error
			if err == nil {
				wrong_trio_counter++
			}
			// Check if the trio is already in the map
			if trio.TBC != "" {
				if _, ok := trio_holder[trio]; ok {
					trio_holder[trio]++
				} else {
					trio_holder[trio] = 1
				}
			}
			//debug.FreeOSMemory()
		}
		if line%100000 == 0 {
			PrintMemUsage()
		}
	}
	return trio_holder

}

func PrintMemUsage() {
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	// For info on each, see: https://golang.org/pkg/runtime/#MemStats
	fmt.Printf("Alloc = %v MiB", bToMb(m.Alloc))
	fmt.Printf("\tTotalAlloc = %v MiB", bToMb(m.TotalAlloc))
	fmt.Printf("\tSys = %v MiB", bToMb(m.Sys))
	fmt.Printf("\tNumGC = %v\n", m.NumGC)
}

func bToMb(b uint64) uint64 {
	return b / 1024 / 1024
}

func main() {
	// Get the file name from the command line
	read1_dir := os.Args[1]
	read2_dir := os.Args[2]
	fileName := PairedRead{read1_dir, read2_dir}
	//PrintMemUsage()
	// Call the reader function
	trio_holder := reader(fileName)
	// Print the files
	for key, value := range trio_holder {
		fmt.Println(key.CellBC, key.Umi, key.TBC, key.RBC, fmt.Sprint(value))
	}
}
