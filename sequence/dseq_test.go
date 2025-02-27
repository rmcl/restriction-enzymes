package sequence

import (
	"reflect"
	"testing"

	"github.com/rmcl/restriction-enzymes/constants"
	"github.com/rmcl/restriction-enzymes/db"
	"github.com/rmcl/restriction-enzymes/enzyme"
)

func TestCutDSeqWatsonStrandOnly(t *testing.T) {
	dSeq1 := NewFromWatsonStrand("NNGGTCTCNCACANNNNCTCT", constants.Linear)

	ez := db.Enzymes["BsaI"]
	fragments := dSeq1.Cut(&ez)

	for _, fragment := range fragments {
		fragment.Print()
	}

	if len(fragments) != 2 {
		t.Errorf("Expected 2 fragments, got %d, Fragments: %v", len(fragments), fragments)
	}

	if fragments[0].Watson != "NNGGTCTCN" {
		t.Errorf("Expected NNGGTCTCN, got %s", fragments[0].Watson)
	}
	if fragments[0].Crick != "NNCCAGAGNGTGT" {
		t.Errorf("Expected NNCCAGAGNGTGT, got %s", fragments[0].Crick)
	}

}

func TestCutDSeqCutSiteOnBothStrands(t *testing.T) {
	dSeq1 := NewFromWatsonStrand("NNGGTCTCNCACANNNNCTCTNGAGACCNN", constants.Linear)

	BsaI := db.Enzymes["BsaI"]
	fragments := dSeq1.Cut(&BsaI)

	for _, fragment := range fragments {
		fragment.Print()
	}

	if len(fragments) != 3 {
		t.Errorf("Expected 3 fragments, got %d, Fragments: %v", len(fragments), fragments)
	}

	expected := []string{
		"NNGGTCTCN", "NNCCAGAGNGTGT", //Frag 1
		"CACANNNN", "NNNNGAGA", // Frag 2
		"CTCTNGAGACCNN", "NCTCTGGNN", // Frag 3
	}

	actual := []string{
		fragments[0].Watson, fragments[0].Crick,
		fragments[1].Watson, fragments[1].Crick,
		fragments[2].Watson, fragments[2].Crick,
	}

	if reflect.DeepEqual(expected, actual) == false {
		t.Errorf("Expected fragments %v, got %v", expected, actual)
	}
}

func TestCutWithRestrictionBatch(t *testing.T) {
	batch := enzyme.NewRestrictionBatch(
		db.Enzymes["BsaI"],
		db.Enzymes["EcoRI"],
	)

	dSeq1 := NewFromWatsonStrand("AATAGACAGAATTCGATTCACCAGAGGTCTCATAGACAAACCAGAGAAAAAAAA", constants.Linear)

	fragments := dSeq1.Cut(&batch)

	actual := []string{
		fragments[0].Watson, fragments[0].Crick,
		fragments[1].Watson, fragments[1].Crick,
		fragments[2].Watson, fragments[2].Crick,
	}

	expected := []string{
		"AATAGACAG", "TTATCTGTCTTAA",
		"AATTCGATTCACCAGAGGTCTCA", "GCTAAGTGGTCTCCAGAGTATCT",
		"TAGACAAACCAGAGAAAAAAAA", "GTTTGGTCTCTTTTTTTT",
	}

	if reflect.DeepEqual(expected, actual) == false {
		t.Errorf("Expected fragments %v, got %v", expected, actual)
	}
}

func TestCutWithRestrictionBatchCircularCutLoop(t *testing.T) {
	// GGTCTC cut site should wrap around
	dSeq1 := NewFromWatsonStrand("AAAGGTCTCNCACANNNNCCAA", constants.Circular)

	// AAAGGTCTCN
	// TTTCCAGAGNGTGT

	// CACANNNNCCAA
	//     NNNNGGTT

	batch := enzyme.NewRestrictionBatch(
		db.Enzymes["BsaI"],
	)

	results := dSeq1.Cut(&batch)

	if len(results) != 2 {
		t.Errorf("Expected 2 fragment, got %d", len(results))
		t.Fail()
	}

	if results[0].Watson != "AAAGGTCTCN" {
		t.Errorf("Expected AAAGGTCTCN, got %s", results[0].Watson)
		t.Fail()
	}

	if results[0].Crick != "TTTCCAGAGNGTGT" {
		t.Errorf("Expected TTTCCAGAGNGTGT, got %s", results[0].Crick)
		t.Fail()
	}

	if results[1].Watson != "CACANNNNCCAA" {
		t.Errorf("Expected CACANNNNCCAA, got %s", results[1].Watson)
		t.Fail()
	}
	if results[1].Crick != "NNNNGGTT" {
		t.Errorf("Expected NNNNGGTT, got %s", results[1].Crick)
		t.Fail()
	}

}
