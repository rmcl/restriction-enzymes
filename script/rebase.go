/*
Package script provides functions for downloading and processing REBASE files.

REBASE is a database of restriction enzymes and their recognition sites.
The database is available for download at ftp://ftp.neb.com/pub/rebase.

The REBASE database is distributed as a set of files, each containing a
different type of information. The files are named according to the following
pattern:

	emboss_e.###	- Enzyme names and recognition sites
	emboss_s.###	- Enzyme suppliers
	emboss_r.###	- Enzyme references
	bairoch.###

where ### represents the REBASE version.
*/
package script

import (
	"bufio"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"strings"

	"github.com/bebop/poly/transform"
	"github.com/rmcl/restriction-enzymes/enzyme"
	"github.com/secsy/goftp"
)

/*

	The REBASE file format header describes the columns as follows:
			# Where:
		# name = name of enzyme
		# pattern = recognition site
		# len = length of pattern
		# ncuts = number of cuts made by enzyme
		#         Zero represents unknown
		# blunt = true if blunt end cut, false if sticky
		# c1 = First 5' cut
		# c2 = First 3' cut
		# c3 = Second 5' cut
		# c4 = Second 3' cut
		#
		# Examples:
		# AAC^TGG -> 6 2 1 3 3 0 0
		# A^ACTGG -> 6 2 0 1 5 0 0
		# AACTGG  -> 6 0 0 0 0 0 0
		# AACTGG(-5/-1) -> 6 2 0 1 5 0 0
		# (8/13)GACNNNNNNTCA(12/7) -> 12 4 0 -9 -14 24 19
		#
		# i.e. cuts are always to the right of the given
		# residue and sequences are always with reference to
		# the 5' strand.
		# Sequences are numbered ... -3 -2 -1 1 2 3 ... with
		# the first residue of the pattern at base number 1.
*/

func convertToInts(intStrings ...string) ([]int, error) {
	intArray := make([]int, len(intStrings))
	for index, val := range intStrings {
		i, err := strconv.Atoi(val)
		if err != nil {
			return intArray, errors.New("error converting string to int")
		}
		intArray[index] = i
	}
	return intArray, nil
}

/* Process the emboss_e.### file */
func processEnzymeFile(enzymeFp io.Reader, enzymes *map[string]enzyme.Enzyme) error {

	// name = name of enzyme
	// pattern = recognition site
	// len = length of pattern
	// ncuts = number of cuts made by enzyme
	// blunt = true if blunt end cut, false if sticky
	// First 5' cut
	// First 3' cut
	// Second 5' cut
	// Second 3' cut
	enzymePattern := regexp.MustCompile(
		`^\s*(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)$`)

	scanner := bufio.NewScanner(enzymeFp)
	for scanner.Scan() {
		line := scanner.Text()

		// Skip comment lines
		if strings.HasPrefix(line, "#") {
			continue
		}

		matches := enzymePattern.FindStringSubmatch(line)
		if len(matches) != 10 {
			continue
		}

		intMatches, err := convertToInts(
			matches[3], // Length
			matches[6], // 5' Cut Position
			matches[7], // 3' Cut Position
			matches[8], // 2nd 5' Cut Position
			matches[9]) // 2nd 3' Cut Position
		if err != nil {
			return err
		}

		// Process Number of Cuts
		var cuts enzyme.EnzymeNumberOfCuts
		switch matches[4] {
		case "0":
			cuts = enzyme.UnknownCuts
		case "1":
			cuts = enzyme.OneCut
		case "2":
			cuts = enzyme.TwoCuts
		default:
			return fmt.Errorf("invalid number of cuts: %s", matches[4])
		}

		// Process Blunt Value
		var cutType enzyme.EnzymeCutType
		switch matches[5] {
		case "0":
			cutType = enzyme.StickyEnd
		case "1":
			cutType = enzyme.BluntEnd
		default:
			return fmt.Errorf("invalid cut type: %s", matches[5])
		}

		enzymeRecord := enzyme.Enzyme{
			Name: matches[1],
			Site: matches[2],

			Length:    intMatches[0],
			RegexpFor: regexp.MustCompile("(?i)" + matches[2]),
			RegexpRev: regexp.MustCompile("(?i)" + transform.ReverseComplement(matches[2])),

			NumberOfCuts: cuts,
			CutType:      cutType,

			Substrate: "DNA",

			FivePrimeCutSite:   intMatches[1],
			ThreePrimeCutSite:  intMatches[2],
			FivePrimeCutSite2:  intMatches[3],
			ThreePrimeCutSite2: intMatches[4],
		}

		if _, ok := (*enzymes)[enzymeRecord.Name]; ok {
			fmt.Printf("duplicate enzyme name: %s", enzymeRecord.Name)

			// For some reason there is at least one duplicate -- HpyUM037X
			// Ignoring and going on with the rest for now.
			continue
		}

		(*enzymes)[enzymeRecord.Name] = enzymeRecord

	}

	return nil
}

/* Process the bairoch.### file

The file is composed of many records separated by "//". Each record is composed of

	Record Example

	ID   AaaI
	ET   R2
	AC   RB00001;
	OS   Acetobacter aceti ss aceti
	PT   XmaIII
	RS   CGGCCG, 1;
	CR   .
	RN   [1]
	RA   Tagami H., Tayama K., Tohyama T., Fukaya M., Okumura H., Kawamura Y.,
	RA   Horinouchi S., Beppu T.;
	RL   FEMS Microbiol. Lett. 56:161-166(1988).
	//
*/

func getNextBairochRecord(scanner *bufio.Scanner) (map[string][]string, error) {
	record := make(map[string][]string)

	lineRegex := regexp.MustCompile(`^\s*(\S+)\s+(.*)$`)

	for scanner.Scan() {
		line := scanner.Text()

		// Skip blanks
		if line == "" {
			continue
		}

		// Skip comment lines
		if strings.HasPrefix(line, "CC") {
			continue
		}

		// Check for end of record
		if line == "//" {
			return record, nil
		}

		matches := lineRegex.FindStringSubmatch(line)
		if len(matches) != 3 {
			return nil, fmt.Errorf("invalid line: %s", line)
		}

		if _, ok := record[matches[1]]; ok {
			record[matches[1]] = append(record[matches[1]], matches[2])
		} else {
			record[matches[1]] = []string{matches[2]}
		}

	}
	return nil, io.EOF
}

func processBairochFile(bairochFp io.Reader, enzymes *map[string]enzyme.Enzyme) error {

	scanner := bufio.NewScanner(bairochFp)
	for {
		record, err := getNextBairochRecord(scanner)
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}

		enzymeId := record["ID"][0]

		if _, ok := (*enzymes)[enzymeId]; !ok {
			fmt.Printf("enzyme present in bairoch file, but not in enzyme definitions: '%s'", enzymeId)

			// For some reason there is at least one enzyme present in the bairoch file
			// but not in the enzyme definitions. Ignoring and going on with the rest for now.
			// Example: AaaI
			continue
		}
		enzyme := (*enzymes)[enzymeId]

		regexpId := regexp.MustCompile(`RB(\d+)`)

		rebaseIdMatch := regexpId.FindStringSubmatch(record["AC"][0])
		if len(rebaseIdMatch) != 2 {
			return fmt.Errorf("invalid REBASE ID: %s: '%s'", record["AC"][0], rebaseIdMatch)
		}
		enzyme.RebaseId, err = strconv.Atoi(rebaseIdMatch[1])
		if err != nil {
			return err
		}
		enzyme.Uri = fmt.Sprintf("https://identifiers.org/rebase:%d", enzyme.RebaseId)
		enzyme.References = record["RA"]

		// Put the modified enzyme back in map
		(*enzymes)[enzymeId] = enzyme
	}

	return nil
}

func ProcessRebaseFiles(
	rebaseInputDir,
	version string,
) (map[string]enzyme.Enzyme, error) {
	enzymes := make(map[string]enzyme.Enzyme)

	enzymeFilePath := filepath.Join(rebaseInputDir, "emboss_e."+version)
	//supplierFilePath := filepath.Join(rebaseDir, "emboss_s."+version)
	//referenceFilePath := filepath.Join(rebaseDir, "emboss_r."+version)
	bairochFilePath := filepath.Join(rebaseInputDir, "bairoch."+version)

	eFP, err := os.Open(enzymeFilePath)
	if err != nil {
		return nil, err
	}

	err = processEnzymeFile(eFP, &enzymes)
	if err != nil {
		return nil, err
	}

	bFP, err := os.Open(bairochFilePath)
	if err != nil {
		return nil, err
	}
	err = processBairochFile(bFP, &enzymes)
	if err != nil {
		return nil, err
	}

	return enzymes, nil
}

func WriteEnzymeJSON(enzymes []enzyme.Enzyme, outputFilePath string) error {
	// Convert enzymes to JSON
	jsonData, err := json.Marshal(enzymes)
	if err != nil {
		return err
	}

	// Make directories if they don't exist
	err = os.MkdirAll(filepath.Dir(outputFilePath), 0755)
	if err != nil {
		return err
	}

	err = os.WriteFile(outputFilePath, jsonData, 0644)
	if err != nil {
		return err
	}

	return nil
}

/* Retrieve Files From FTP */

var rebaseFTPHost = "ftp.neb.com"

func rebaseFtpPaths(version string) []string {
	paths := []string{
		"/pub/rebase/emboss_e.###",
		"/pub/rebase/emboss_s.###",
		"/pub/rebase/emboss_r.###",
		"/pub/rebase/bairoch.###",
	}

	for i, path := range paths {
		paths[i] = strings.Replace(path, "###", version, 1)
	}

	return paths
}

func retrieveFileFromFTP(ftpUrl string, downloadDir string) error {
	fileName := filepath.Base(ftpUrl)
	outputFilePath := filepath.Join(downloadDir, fileName)

	fmt.Println("Downloading file from FTP: ", fileName)

	ftpConnection, err := goftp.Dial(rebaseFTPHost)
	if err != nil {
		return err
	}
	defer ftpConnection.Close()

	// Create a local file to save to
	file, err := os.Create(outputFilePath)
	if err != nil {
		return err
	}
	defer file.Close()

	// Get the file from the FTP server
	err = ftpConnection.Retrieve(ftpUrl, file)
	if err != nil {
		return err
	}

	return nil
}

func RetrieveRebaseFiles(version string, downloadDir string) error {
	paths := rebaseFtpPaths(version)

	for _, path := range paths {
		err := retrieveFileFromFTP(path, downloadDir)
		if err != nil {
			return err
		}
	}
	return nil
}
