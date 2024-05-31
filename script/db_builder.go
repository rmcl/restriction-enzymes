package script

import (
	"bytes"
	"os"
	"strings"
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

const enzymeRecordTemplate = `"{{.Enzyme.Name}}":{` +
	`Name:"{{.Enzyme.Name}}",` +
	`Site:"{{.Enzyme.Site}}",` +
	`Length:{{.Enzyme.Length}},` +
	`Substrate:"{{.Enzyme.Substrate}}",` +
	`RegexpFor:regexp.MustCompile("{{.Enzyme.RegexpFor}}"),` +
	`RegexpRev:regexp.MustCompile("{{.Enzyme.RegexpRev}}"),` +
	`OverhangLength:{{.Enzyme.OverhangLength}},` +
	`OverhangSequence:"{{.Enzyme.OverhangSequence}}",` +
	`NumberOfCuts:{{.Enzyme.NumberOfCuts}},` +
	`CutType:"{{.Enzyme.CutType}}",` +
	`FivePrimeCutSite:{{.Enzyme.FivePrimeCutSite}},` +
	`ThreePrimeCutSite:{{.Enzyme.ThreePrimeCutSite}},` +
	`FivePrimeCutSite2:{{.Enzyme.FivePrimeCutSite2}},` +
	`ThreePrimeCutSite2:{{.Enzyme.ThreePrimeCutSite2}},` +
	`RebaseId:{{.Enzyme.RebaseId}},` +
	`InactivationTemperature:{{.Enzyme.InactivationTemperature}},` +
	`OptimalTemperature:{{.Enzyme.OptimalTemperature}},` +
	`Uri:"{{.Enzyme.Uri}}",` +
	`Suppliers: []string{ {{.Suppliers}} },` +
	`References: []string{ {{.References}} },` +
	`},`

type SerializableEnzymeRecord struct {
	Enzyme     enzyme.Enzyme
	References string
	Suppliers  string
}

func NewSerializableEnzymeRecord(e enzyme.Enzyme) SerializableEnzymeRecord {
	references := ""
	if len(e.References) > 0 {
		references = "\"" + strings.Join(e.References, "\", \"") + "\""
	}
	suppliers := ""
	if len(e.Suppliers) > 0 {
		suppliers = "\"" + strings.Join(e.Suppliers, "\", \"") + "\""
	}

	return SerializableEnzymeRecord{
		Enzyme:     e,
		References: references,
		Suppliers:  suppliers,
	}
}

type DBFileTemplateData struct {
	Enzymes string
}

func serializeEnzymeRecord(enzyme enzyme.Enzyme) (string, error) {
	recordTemplate, err := template.New("source").Parse(enzymeRecordTemplate)
	if err != nil {
		return "", err
	}

	var output bytes.Buffer

	err = recordTemplate.Execute(
		&output,
		NewSerializableEnzymeRecord(enzyme))

	if err != nil {
		return "", err
	}

	return output.String(), nil
}

func createEnzymeFileString(enzymes map[string]enzyme.Enzyme) (string, error) {
	records := ""
	for _, enzyme := range enzymes {
		recordString, err := serializeEnzymeRecord(enzyme)
		if err != nil {
			return "", err
		}

		records += "\t" + recordString + "\n"
	}

	fileTemplate, err := template.New("source").Parse(dbFileTemplate)
	if err != nil {
		return "", err
	}

	var buffer bytes.Buffer
	fileTemplate.Execute(&buffer, DBFileTemplateData{
		Enzymes: records,
	})

	return buffer.String(), nil
}

func CreateGoEnzymeDBFile(enzymes map[string]enzyme.Enzyme, outputFilePath string) error {
	fileString, err := createEnzymeFileString(enzymes)
	if err != nil {
		return err
	}

	err = os.WriteFile(outputFilePath, []byte(fileString), 0644)
	if err != nil {
		return err
	}

	return nil
}
