package main

import (
	"fmt"
	"os"

	"github.com/rmcl/restriction-enzymes/script"
)

func getTempDir() string {
	tempDir := os.TempDir()
	return tempDir
}

func main() {

	version := "405"

	tempDir := getTempDir()

	fmt.Println("Downloading files from REBASE FTP to: ", tempDir)

	err := script.RetrieveRebaseFiles(version, tempDir)
	if err != nil {
		panic(err)
	}

	data, err := script.ProcessRebaseFiles(tempDir, version)
	if err != nil {
		panic(err)
	}

	outputGoFilePath := "./db/enzyme_db.go"
	outputGoSupplierFilePath := "./db/suppliers.go"
	err = script.CreateGoEnzymeDBFile(data, outputGoFilePath)
	if err != nil {
		panic(err)
	}
	err = script.CreateGoEnzymeSupplierFile(data, outputGoSupplierFilePath)
	if err != nil {
		panic(err)
	}
}
