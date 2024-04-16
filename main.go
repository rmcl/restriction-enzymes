package main

import (
	"os"

	"github.com/rmcl/restriction-enzymes/script"
)

func getTempDir() string {
	tempDir := os.TempDir()
	return tempDir
}

func main() {

	version := "404"

	tempDir := getTempDir()
	err := script.RetrieveRebaseFiles(version, tempDir)
	if err != nil {
		panic(err)
	}

	enzymes, err := script.ProcessRebaseFiles(tempDir, version)
	if err != nil {
		panic(err)
	}

	outputGoFilePath := "./db/enzyme_db.go"
	script.CreateGoEnzymeDBFile(enzymes, outputGoFilePath)
}
