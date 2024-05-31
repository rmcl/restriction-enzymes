package script

import (
	"bytes"
	"os"
	"sort"
	"text/template"

	"github.com/rmcl/restriction-enzymes/enzyme"
)

const dbFileTemplate = `
/*
Package db contains a map of enzyme names to enzyme.Enzyme structs pulled from the REBASE database.

THIS FILE IS AUTO-GENERATED. DO NOT MODIFY THIS FILE MANUALLY.

To update the database, run the script in the script directory.

*/
package db

import (
	"regexp"
	"github.com/rmcl/restriction-enzymes/enzyme"
)

var Enzymes = map[string]enzyme.Enzyme{
{{.Enzymes}}
}
`

const supplierFileTemplate = `
/*
This file contains a list of Reagenet Suppliers and the restriction enzymes they supply. The data comes from the REBASE database.

THIS FILE IS AUTO-GENERATED. DO NOT MODIFY THIS FILE MANUALLY.

To update the database, run the script in the script directory.

*/
package db

type SupplierRecord struct {
	Id   string
	Name string
	Enzymes      []string
}

var Suppliers = []SupplierRecord{
{{.Suppliers}}
}
`
const supplierRecordTemplate = `{` +
	`Id:"{{.SupplierId}}", ` +
	`Name:"{{.SupplierName}}", ` +
	`Enzymes: []string{ {{formatStringList .Enzymes}} },` +
	`},`

const enzymeRecordTemplate = `"{{.Name}}":{` +
	`Name:"{{.Name}}",` +
	`Site:"{{.Site}}",` +
	`Length:{{.Length}},` +
	`Substrate:"{{.Substrate}}",` +
	`RegexpFor:regexp.MustCompile("{{.RegexpFor}}"),` +
	`RegexpRev:regexp.MustCompile("{{.RegexpRev}}"),` +
	`OverhangLength:{{.OverhangLength}},` +
	`OverhangSequence:"{{.OverhangSequence}}",` +
	`NumberOfCuts:{{.NumberOfCuts}},` +
	`CutType:"{{.CutType}}",` +
	`FivePrimeCutSite:{{.FivePrimeCutSite}},` +
	`ThreePrimeCutSite:{{.ThreePrimeCutSite}},` +
	`FivePrimeCutSite2:{{.FivePrimeCutSite2}},` +
	`ThreePrimeCutSite2:{{.ThreePrimeCutSite2}},` +
	`RebaseId:{{.RebaseId}},` +
	`InactivationTemperature:{{.InactivationTemperature}},` +
	`OptimalTemperature:{{.OptimalTemperature}},` +
	`Uri:"{{.Uri}}",` +
	`References: []string{ {{formatStringList .References}} },` +
	`},`

// formatStringList takes a list of strings and returns a string with each string in the list
// wrapped in double quotes and separated by commas.
func formatStringList(stringList []string) string {
	output := ""
	for _, s := range stringList {
		output += "\"" + s + "\", "
	}

	if len(output) > 2 {
		output = output[:len(output)-2]
	}
	return output
}

type DBFileTemplateData struct {
	Enzymes string
}

type SupplierFileTemplateData struct {
	Suppliers string
}

// Give a template and data, parse the template and return the result as a string.
func parseTemplate(templateStr string, data interface{}) (string, error) {
	funcMap := template.FuncMap{
		"formatStringList": formatStringList,
	}

	recordTemplate, err := template.New("source").Funcs(funcMap).Parse(templateStr)
	if err != nil {
		return "", err
	}

	var output bytes.Buffer

	err = recordTemplate.Execute(
		&output,
		data)

	if err != nil {
		return "", err
	}

	return output.String(), nil
}

func createEnzymeFile(enzymes map[string]enzyme.Enzyme) (string, error) {
	records := ""
	for _, enzyme := range enzymes {
		recordString, err := parseTemplate(enzymeRecordTemplate, enzyme)
		if err != nil {
			return "", err
		}

		records += "\t" + recordString + "\n"
	}

	return parseTemplate(dbFileTemplate, DBFileTemplateData{
		Enzymes: records,
	})
}

type SupplierRecord struct {
	SupplierId   string
	SupplierName string
	Enzymes      []string
}

func buildEnzymeSupplierRecords(rebaseData *RebaseData) ([]SupplierRecord, error) {

	// Create a mapping of Supplier 1 letter code to a list of enzyme names
	enzymeNamesBySupplierId := map[string][]string{}
	for _, reference := range rebaseData.References {
		for _, supplierId := range reference.Suppliers {
			enzymeNamesBySupplierId[supplierId] = append(
				enzymeNamesBySupplierId[supplierId],
				reference.EnzymeName)
		}
	}

	supplierRecords := []SupplierRecord{}

	for supplierId, enzymeNames := range enzymeNamesBySupplierId {
		sr := SupplierRecord{
			SupplierId:   supplierId,
			SupplierName: rebaseData.Suppliers[supplierId],
			Enzymes:      enzymeNames,
		}
		supplierRecords = append(supplierRecords, sr)
	}

	// Sort the supplier records by supplier name
	sort.Slice(supplierRecords, func(i, j int) bool {
		return supplierRecords[i].SupplierId < supplierRecords[j].SupplierId
	})

	return supplierRecords, nil
}

func createEnzymeSupplierFile(rebaseData *RebaseData) (string, error) {
	supplierRecords, err := buildEnzymeSupplierRecords(rebaseData)
	if err != nil {
		return "", err
	}

	records := ""
	for _, supplierRecord := range supplierRecords {
		recordString, err := parseTemplate(supplierRecordTemplate, supplierRecord)
		if err != nil {
			return "", err
		}

		records += "\t" + recordString + "\n"
	}

	return parseTemplate(supplierFileTemplate, SupplierFileTemplateData{
		Suppliers: records,
	})
}

func CreateGoEnzymeDBFile(rebaseData *RebaseData, outputFilePath string) error {
	fileString, err := createEnzymeFile(rebaseData.Enzymes)
	if err != nil {
		return err
	}

	err = os.WriteFile(outputFilePath, []byte(fileString), 0644)
	if err != nil {
		return err
	}

	return nil
}

func CreateGoEnzymeSupplierFile(rebaseData *RebaseData, outputFilePath string) error {
	fileString, err := createEnzymeSupplierFile(rebaseData)
	if err != nil {
		return err
	}

	err = os.WriteFile(outputFilePath, []byte(fileString), 0644)
	if err != nil {
		return err
	}

	return nil
}
