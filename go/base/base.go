package base

//import "fmt"

// TODO

type Atom struct {
	name string
	x    float64
	y    float64
	z    float64
}

type XYZ struct {
	file          string
	natom         int
	nspecies      int
	atoms         []Atom
	specie_labels map[string]int
}

func (a Atom) SetName(name string) {
	a.name = name
}

func (a Atom) Setx(x float64) {
	a.x = x
}

func (a Atom) Sety(y float64) {
	a.y = y
}

func (a Atom) Setz(z float64) {
	a.z = z
}

func (a XYZ) GetInfo() {
	//
}

func (a XYZ) SetSpeciesNumber() {
	//
}

func (a XYZ) Tofdf(fname string) {
	//
}

func (a XYZ) Update(newfile string) {
	//
}
