package main

import (
	"fmt"
	"math"
	"os"

	"github.com/veandco/go-sdl2/sdl"
)

var names = []string{
	"sun",
	"earth",
	"ariel",
	"callisto",
	"moon",
	"deimos",
	"dione",
	"enceladus",
	"europa",
	"ganymede",
	"hyperion",
	"iapetus",
	"io",
	"jupiter",
	"mars",
	"mercury",
	"mimas",
	"miranda",
	"neptune",
	"nereid",
	"oberon",
	"phobos",
	"phoebe",
	"pluto",
	"rhea",
	"saturn",
	"tethys",
	"titan",
	"triton",
	"umbriel",
	"uranus",
	"venus",
}

type bri int

const (
	br0 bri = iota
	br1
	br2
	br3
)

var points = []bri{
	br3,
	br2,
	br0,
	br2,
	br1,
	br0,
	br1,
	br0,
	br1,
	br2,
	br0,
	br1,
	br1,
	br3,
	br2,
	br1,
	br0,
	br0,
	br3,
	br0,
	br1,
	br0,
	br0,
	br1,
	br1,
	br3,
	br1,
	br2,
	br1,
	br0,
	br3,
	br2,
}

var ppar = []int{
	0, 0, 036, 015, 01, 016, 031, 031,
	015, 015, 031, 031, 015, 0, 0, 0,
	031, 036, 0, 022, 036, 016, 031, 0,
	031, 0, 031, 031, 022, 036, 0, 0,
}

// format for hand-rolled floating point used in PDP-7 version
type flt struct {
	exp int8 // exponent
	m1  int  // signed 36-bit mantissa (0-18)
	m2  int  // 18-36
}

var prsq [32]float64
var prsq_ = [...]flt{
	{016, 0272245, 075341},
	{1, 0200000, 0},
	{-07, 0244122, 0506362},
	{-02, 0251477, 0620663},
	{-03, 0230761, 0127762},
	{-025, 0320300, 054474},
	{-06, 0324134, 0124211},
	{-010, 0335416, 0541570},
	{-04, 0371372, 0},
	{-06, 0247430, 0},
	{-011, 0311150, 0},
	{-05, 0302622, 0},
	{-02, 0256475, 0},
	{07, 0376733, 0},
	{-01, 0221530, 0},
	{-02, 0235142, 0},
	{-010, 0217266, 0},
	{-011, 0274361, 0},
	{04, 0365471, 0},
	{-012, 0227176, 0},
	{-06, 0342454, 0},
	{-023, 0326340, 0},
	{-013, 0326265, 0},
	{-02, 0323774, 0},
	{-05, 0255140, 0},
	{07, 0263573, 0},
	{-06, 0223174, 0},
	{-02, 0251477, 0},
	{-02, 0235142, 0},
	{-05, 0223060, 0},
	{05, 0206115, 0},
	{0, 0362406, 0},
}

var accl [32]float64
var accl_ = [...]flt{
	{0, 0204365, 0},
	{-023, 0320324, 0},
	{-036, 0227207, 0},
	{-017, 0340500, 0},
	{-030, 0210041, 0},
	{-063, 0341666, 0},
	{-034, 0235122, 0},
	{-037, 0247531, 0},
	{-031, 0310316, 0},
	{-027, 0334427, 0},
	{-041, 0315203, 0},
	{-033, 0303403, 0},
	{-030, 0245752, 0},
	{-017, 0201414, 0},
	{-026, 0263753, 0},
	{-026, 0205241, 0},
	{-040, 0256464, 0},
	{-041, 0272051, 0},
	{-017, 0340566, 0},
	{-043, 0275073, 0},
	{-034, 0255345, 0},
	{-060, 0341650, 0},
	{-044, 0341552, 0},
	{-020, 0307762, 0},
	{-033, 0243712, 0},
	{-014, 0233053, 0},
	{-035, 0265543, 0},
	{-027, 0340500, 0},
	{-027, 0210344, 0},
	{-037, 0210777, 0},
	{-017, 0275653, 0},
	{-023, 0252667, 0},
}

var px [32]float64
var px_ = [...]flt{
	{0, 000000, 0},
	{015, 0620356, 0},
	{005, 0360005, 0},
	{010, 0666214, 0},
	{005, 0704053, 0},
	{002, 0347600, 0},
	{005, 0310506, 0},
	{006, 0220622, 0},
	{007, 0310473, 0},
	{004, 0370065, 0},
	{006, 0304101, 0},
	{011, 0676631, 0},
	{006, 0653020, 0},
	{020, 0317202, 0},
	{017, 0644356, 0},
	{013, 0206414, 0},
	{005, 0245346, 0},
	{004, 0222264, 0},
	{023, 0261234, 0},
	{001, 0372225, 0},
	{007, 0646102, 0},
	{001, 0212446, 0},
	{013, 0773152, 0},
	{024, 0274557, 0},
	{004, 0227474, 0},
	{020, 0263122, 0},
	{004, 0333254, 0},
	{010, 0216672, 0},
	{006, 0231142, 0},
	{006, 0212701, 0},
	{023, 0650051, 0},
	{016, 0233751, 0},
}

var py [32]float64
var py_ = [...]flt{
	{000, 0000000, 0},
	{017, 0664054, 0},
	{002, 0662035, 0},
	{010, 0350757, 0},
	{006, 0334771, 0},
	{000, 0267000, 0},
	{006, 0726770, 0},
	{004, 0225752, 0},
	{006, 0201346, 0},
	{010, 0247536, 0},
	{010, 0343277, 0},
	{011, 0314411, 0},
	{006, 0712237, 0},
	{021, 0334656, 0},
	{017, 0342324, 0},
	{016, 0615151, 0},
	{005, 0644257, 0},
	{005, 0622456, 0},
	{024, 0224063, 0},
	{012, 0331314, 0},
	{006, 0640034, 0},
	{001, 0200024, 0},
	{011, 0243161, 0},
	{023, 0772355, 0},
	{007, 0644432, 0},
	{022, 0726324, 0},
	{006, 0260740, 0},
	{010, 0600213, 0},
	{006, 0237476, 0},
	{005, 0674734, 0},
	{023, 0616334, 0},
	{016, 0327155, 0},
}

var pw [32]float64

var pww [32]float64
var pww_ = [...]flt{
	{0000, 0000000, 0000000},
	{-054, 0663265, 0376074},
	{-036, 0743326, 0460356},
	{-043, 0647730, 0444215},
	{-045, 0767246, 0341205},
	{-034, 0745027, 0221674},
	{-036, 0702670, 0530661},
	{-034, 0702340, 0273047},
	{-037, 0747771, 0100452},
	{-041, 0743411, 0732756},
	{-044, 0716603, 0200021},
	{-050, 0755641, 0517072},
	{-035, 0751602, 0332677},
	{-072, 0631066, 0300145},
	{-056, 0712656, 0015171},
	{-050, 0701167, 0507203},
	{-033, 0715124, 0202507},
	{-034, 0665024, 0365605},
	{-073, 0730064, 0776551},
	{-054, 0667362, 0431776},
	{-043, 0775760, 0312631},
	{-030, 0740613, 0034530},
	{-055, 0635657, 0370276},
	{-074, 0677234, 0315321},
	{-037, 0617102, 0655555},
	{-066, 0723370, 0773672},
	{-035, 0714526, 0724272},
	{-043, 0667633, 0402706},
	{-040, 0636053, 0440472},
	{-037, 0650025, 0224325},
	{-071, 0717663, 0701773},
	{-052, 0754612, 0304722},
}

func flt2float(flt flt) float64 {
	a := flt.exp
	b := flt.m1
	c := flt.m2
	s := 1 - 2*(b>>17)
	b &= 0377777
	m := 0.0
	i := 1
	for b != 0 {
		m += float64(b>>16) / float64(int(1)<<i)
		b <<= 1
		b &= 0377777
		i++
	}
	i = 18
	for c != 0 {
		m += float64(c>>17) / float64(int(1)<<i)
		c <<= 1
		c &= 0777777
		i++
	}
	if a > 0 {
		return float64(s) * m * float64(int(1)<<a)
	}
	return float64(s) * m / float64(int(1)<<-a)
}

var (
	rpar, dpar   float64
	ax, ay, maxa float64
	maxj         int

	sphi, cphi float64
	absx, absy float64
	spx, spy   float64

	shipx, shipy, x, y, ox, oy    float64
	lanflg, goflg, forflg, bacflg bool

	horizv float64

	inflg, grvflg bool
	crflg         bool

	par    int
	stheta = 1.0
	ctheta = 0.0
	scale  = 0
)

// Constants
const (
	nplan  = 32
	vscale = 6
	ascale = -0.5
	fardst = 32768
	thrs   = 2
	accflg = false
)

var (
	locflg       bool
	locpar       int
	crash        float64
	sdphi, cdphi float64

	xpos = 0
	ypos = 0
	sz   = 1
	br   = 1
	blnk = false
	show = true
)

const screenWidth = 1024
const screenHeight = 1024

var renderer *sdl.Renderer

func dsetx(x int) {
	xpos = x
}

func dsety(y int) {
	ypos = y
}

func dscale(s int) {
	sz = 1 << s
}

func intens(i bri) {
	br = int(i) + 1
}

func blink(s bool) {
	blnk = s
}

func vec(x, y int, vis bool) {
	if vis {
		renderer.SetDrawColor(uint8(255*br/4), uint8(255*br/4), uint8(255*br/4), 255)
		renderer.DrawLine(int32(xpos), int32(screenHeight-ypos), int32(xpos+sz*x), int32(screenHeight-(ypos+sz*y)))
	}
	xpos += sz * x
	ypos += sz * y
}

func vecx(x int) {
	vec(x, 0, true)
}

func vecy(y int) {
	vec(0, y, true)
}

type dir int

const (
	N dir = iota
	NE
	E
	SE
	S
	SW
	W
	NW
)

func incr(n int, dir dir) {
	xtab := [...]int{0, 1, 1, 1, 0, -1, -1, -1}
	ytab := [...]int{1, 1, 0, -1, -1, -1, 0, 1}
	dx := xtab[dir]
	dy := ytab[dir]
	for i := 0; i < n; i++ {
		vec(dx, dy, true)
	}
}

func dchar(c byte) {
	vis := false
	ch := font[c-32]
	x := 0
	y := 3
	for _, ij := range ch {
		if ij == 0xff {
			break
		}
		if ij == 0x0f {
			vis = false
		} else if ij == 0x1f {
			vis = true
		} else {
			x2 := ij & 0xf
			y2 := 0xf - (ij >> 4)
			dx := int(x2) - x
			dy := int(y2) - y
			vec(dx, dy, vis)
			x += dx
			y += dy
		}
	}
	vec(13-x, 3-y, false)
}

func chars(s string) {
	if blnk && !show {
		return
	}
	for _, r := range s {
		if r == 0 {
			break
		}
		dchar(byte(r))
	}
}

func inscr(a int) int {
	a = 383 - a
	if a < 0 {
		return 0
	}
	a -= 768
	if a >= 0 {
		return 0
	}
	return a
}

func rotx() {
	s := absx*stheta - absy*ctheta
	if scale > 0 {
		spx = s / float64(int(1)<<scale)
	} else {
		spx = s * float64(int(1)<<-scale)
	}
}

var sca [4]byte

func dssca() {
	chars(string(sca[:]))
	dscale(0)
	dsetx(127)
	dsety(250)
	vecx(768)
	dsetx(895)
	dsety(255)
	vecy(768)
	dsetx(895)
	dsety(1023)
	vecx(-768)
	dsetx(127)
	dsety(1023)
	vecy(-768)
	dsetx(127)
	dsety(255)
	vecx(768)
	dsetx(511)
	dsety(255)
	vecy(767)
	dsetx(127)
	dsety(639)
	vecx(767)
}

func dspsca() {
	var s int
	if scale < 0 {
		s = -scale
		sca[0] = '-'
	} else {
		s = scale
		sca[0] = '+'
	}
	sca[1] = '0' + byte(s/10)
	sca[2] = '0' + byte(s%10)
}

func rotate(dir bool) {
	ftmp1 := cphi * sdphi
	if !dir {
		ftmp1 = -ftmp1
	}
	ftmp2 := sphi*cdphi + ftmp1
	ftmp1 = sphi * sdphi
	if dir {
		ftmp1 = -ftmp1
	}
	cphi = cphi*cdphi + ftmp1
	sphi = ftmp2
}

var quit = false
var pbson []uint8

func contrl(e sdl.Event) {
	if e != nil {
		switch e.GetType() {
		case sdl.KEYDOWN:
			ke, _ := e.(*sdl.KeyboardEvent)
			switch ke.Keysym.Sym {
			case sdl.K_1:
				quit = true
			case sdl.K_2:
				goflg = false
				crflg = false
			case sdl.K_7, sdl.K_DOWN:
				scale++
				dspsca()
			case sdl.K_8, sdl.K_UP:
				scale--
				dspsca()
			}
		}
	}
	forflg = false
	bacflg = false
	if pbson[sdl.SCANCODE_3] != 0 {
		forflg = true
		if !goflg {
			lanflg = false
		}
	}
	if pbson[sdl.SCANCODE_4] != 0 {
		bacflg = true
		if !goflg {
			lanflg = false
		}
	}
	if pbson[sdl.SCANCODE_5] != 0 || pbson[sdl.SCANCODE_RIGHT] != 0 {
		rotate(false)
	}
	if pbson[sdl.SCANCODE_6] != 0 || pbson[sdl.SCANCODE_LEFT] != 0 {
		rotate(true)
	}
}

func dsplanet(p int) {
	intens(points[p])
	incr(1, NE)
	incr(2, W)
	incr(2, S)
	incr(2, E)
	incr(1, N)
	incr(1, W)
}

var cl = " "

func dispcl() {
	chars(cl)
}

func namedsp() {
	chars(names[par])
}

func displist() {
	dscale(1)
	intens(3)
	blink(true)
	dsetx(800)
	dsety(20)
	dispcl()
	intens(0)
	blink(false)
	dsetx(0)
	dsety(20)
	namedsp()
	dsetx(400)
	dsety(20)
	dssca()
}

// Update planet position
func updpln(p int) {
	// rotate by a (pw[p] = cos(a), pww[p] = sin(a))
	ftmp1 := px[p]
	px[p] = px[p]*pw[p] - py[p]*pww[p]
	py[p] = ftmp1*pww[p] + py[p]*pw[p]
}

func invert(p int) {
	pww[p] = -pww[p]
	updpln(p)
}

func absv(p int) {
	absx = 0
	absy = 0

	absi := p
	for absi != 0 {
		invert(absi)
		absx -= px[absi]
		absy -= py[absi]

		invert(absi)
		absx += px[absi]
		absy += py[absi]

		absi = ppar[absi]
	}
}

// Get absolute planet position
func absxy(p int) {
	// get distance from planet to sun
	absi := p
	absx = 0
	absy = 0
	for absi != 0 {
		absx += px[absi]
		absy += py[absi]
		absi = ppar[absi]
	}
}

// Set absolute ship position
func shipxy() {
	shipx = -(absx + x)
	shipy = -(absy + y)
}

// Update ship acceleration from planet gravity
func updacc(p int) {
	// set absx and absy to distance from ship to planet
	if p != par {
		absxy(p)
		absx += shipx
		absy += shipy
	} else {
		absx = -x
		absy = -y
	}

	dtmp1 := absy*absy + absx*absx

	dpar = math.Sqrt(dtmp1)

	if p == par {
		ftmp1 := x - ox
		ftmp2 := oy - y
		horizv = (ftmp2*x + ftmp1*y) / dpar
		// if distance from planet less than radius
		if dpar < rpar {
			if !lanflg {
				if ftmp1*ftmp1+ftmp2*ftmp2 > crash {
					crflg = true
					goflg = true
				}
				lanflg = true
				ftmp1 = rpar / dpar
				// reset x and y to the edge of the planet
				x *= ftmp1
				ox = x
				y *= ftmp1
				oy = y
				absxy(par)
				shipxy()
				updacc(par)
			}
		}
	}

	/* upda5 */
	// if the planet is too far, set a flag to not draw its border
	// and ignore its gravity
	if dpar > fardst {
		// unless it's the sun
		if p != 0 {
			grvflg = true
			return
		}
	}
	grvflg = false

	// a = GM / r^2
	a := accl[p] / dtmp1
	if a > maxa {
		maxa = a
		maxj = p
	}

	// set the new acceleration
	ftmp1 := a / dpar
	ax += ftmp1 * absx
	ay += ftmp1 * absy
}

// Draw circle
func surf(nt, setx, sety int, wx, wy, v, vv float64) bool {
	tsetx := -setx
	dsetx(895 + setx)
	tsety := -sety
	dsety(1023 + sety)
	twx := wx
	twy := wy
	v = -v
	for i := 0; i < nt; i++ {
		// rotate by a (vv = cos(a), v = sin(a))
		ftmp2 := vv*twx - v*twy
		twy = vv*twy + v*twx
		res := inscr(int(twy + spy))
		if res == 0 {
			return false
		}
		dely := res + tsety
		tsety -= dely
		twx = ftmp2
		res = inscr(int(twx + spx))
		if res == 0 {
			return false
		}
		delx := res + tsetx
		tsetx -= delx
		vec(delx, dely, true)
	}
	return true
}

// Draw planet circle
func drcirc(p int) {
	if grvflg {
		return
	}
	pr := math.Sqrt(prsq[p])
	var dtmp1 float64
	if scale > 0 {
		dtmp1 = pr / float64(int(1)<<scale)
	} else {
		dtmp1 = pr * float64(int(1)<<-scale)
	}

	if dtmp1 < thrs {
		return
	}

	if scale > 0 {
		dpar = dpar / float64(int(1)<<scale)
	} else {
		dpar = dpar * float64(int(1)<<-scale)
	}
	dtmp2 := (dpar - dtmp1) / dpar

	wy := dtmp2 * spy
	res := inscr(int(wy))
	if res == 0 {
		return
	}
	sety := res

	if !inflg {
		rotx()
	}

	wx := dtmp2 * spx
	res = inscr(int(wx))
	if res == 0 {
		return
	}
	setx := res

	wy -= spy
	wx -= spx

	var narcs int
	dtmp1 = dtmp1 * math.Pi / 10
	if dtmp1 > 400 {
		narcs = 400
	} else {
		narcs = int(dtmp1)
		if narcs < 20 {
			narcs = 0
		}
		narcs += 20
		dtmp1 = float64(narcs)
	}

	// v ~= sin(x) ~= x (at x ~ 0)
	// vv ~= cos(x) = sqrt(1 - sin2(x)) ~= 1 - sin2(x)/2
	// (taylor series expansion sqrt(1 - x^2) = 1 - x^2/2 + ...)
	v := (2 * math.Pi) / dtmp1
	vv := 1 - (v*v)/2

	intens(0)
	if surf(narcs, setx, sety, wx, wy, v, vv) {
		return
	}
	surf(narcs, setx, sety, wx, wy, -v, vv)
}

// Display planet
func displa(p int) {
	if p == locpar {
		if locflg {
			stheta = (sphi*absx + cphi*absy) / dpar
			ctheta = (cphi*absx - sphi*absy) / dpar
		} else {
			stheta = sphi
			ctheta = cphi
		}
	}
	f := absy*stheta + absx*ctheta
	if scale > 0 {
		spy = f / float64(int(1)<<scale)
	} else {
		spy = f * float64(int(1)<<-scale)
	}
	inflg = false
	res := inscr(int(spy))
	if res != 0 {
		dsety(1023 + res)
		rotx()
		inflg = true
		res = inscr(int(spx))
		if res != 0 {
			dsetx(895 + res)
			dsplanet(p)
		}
	}
	drcirc(p)
}

// Update ship thrust acceleration and position
func updshp() {
	if forflg || bacflg {
		var a float64
		if forflg && bacflg {
			a = 0
		} else if scale > 0 {
			a = ascale * float64(int(1)<<scale)
		} else {
			a = ascale / float64(int(1)<<-scale)
		}
		if bacflg {
			a = -a
		}
		if forflg && accflg {
			a += maxa
		}
		a = -a
		ax += a * ctheta
		ay += a * stheta
	}

	// move ship by current thrust
	ftmp1 := x*2 - ox + ax
	ox = x
	x = ftmp1
	ftmp1 = y*2 - oy + ay
	oy = y
	y = ftmp1

	// return if we're already set to the planet with strongest gravity
	if par == maxj {
		return
	}

	// update ship position
	absxy(par)
	shipxy()

	// ??
	absv(par)
	ox = absx + x - ox
	oy = absy + y - oy

	// set the current planet to the one with strongest gravity
	par = maxj

	// ??
	absv(par)
	ox = absx - ox
	oy = absy - oy

	// set new relative position
	absxy(par)
	x = -(absx + shipx)
	ox += x
	y = -(absy + shipy)
	oy += y

	// update current planet radius
	rpar = math.Sqrt(prsq[par])

	dspsca()
}

func loop() {
	// set ship position
	absxy(par)
	shipxy()
	// if not game over
	if !goflg {
		ax = 0
		ay = 0
		maxa = 0
		/* loop1 */
		// update acceleration and draw planets
		for i := 0; i < nplan; i++ {
			updacc(i)
			displa(i)
			if i != 0 {
				updpln(i)
			}
		}
		/* loop2 */
		// update ship if we haven't landed
		if !lanflg {
			updshp()
		}
	}
	/* loop3 */
	mul := vscale - scale - 1
	var h float64
	if mul > 0 {
		h = horizv * float64(int(1)<<mul)
	} else {
		h = horizv / float64(int(1)<<-mul)
	}
	res := inscr(int(h))
	if res != 0 {
		dsetx(895 + res)
		dsety(250)
		dsplanet(0)
	}
	/* loop4 */
	if crflg {
		cl = "CL"
	} else if lanflg {
		cl = "L"
	} else {
		cl = " "
	}
}

func main() {
	crash = flt2float(flt{-027, 0200000, 0})
	sdphi = math.Sin(1.2 * math.Pi / 180)
	cdphi = math.Cos(1.2 * math.Pi / 180)

	for i := 0; i < nplan; i++ {
		prsq[i] = flt2float(prsq_[i])
		accl[i] = flt2float(accl_[i])
		px[i] = flt2float(px_[i])
		py[i] = flt2float(py_[i])
		pww[i] = flt2float(pww_[i])
		pw[i] = 1 - (pww[i]*pww[i])/2
	}

	lanflg = true
	crflg = false
	goflg = false
	forflg = false
	bacflg = false
	locflg = false
	locpar = 0

	// start on planet earth
	par = 1

	// set initial position and orientation
	oy = math.Sqrt(prsq[par])
	y = oy
	rpar = oy
	stheta = 1
	sphi = 1
	ctheta = 0
	cphi = 0
	ox = 0
	x = 0

	dspsca()

	if err := sdl.Init(sdl.INIT_VIDEO); err != nil {
		fmt.Fprintf(os.Stderr, "Failed to initialize SDL: %s\n", sdl.GetError())
		os.Exit(1)
	}

	var window *sdl.Window

	var err error
	window, renderer, err = sdl.CreateWindowAndRenderer(screenWidth, screenHeight, sdl.WINDOW_SHOWN)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Failed to create window: %s\n", sdl.GetError())
		sdl.Quit()
		os.Exit(1)
	}

	window.SetTitle("Space Travel")

	pbson = sdl.GetKeyboardState()

	var e sdl.Event
	secs := uint64(0)
	for !quit {
		for {
			e = sdl.PollEvent()
			if e == nil {
				break
			}
			if e.GetType() == sdl.QUIT {
				quit = true
			} else {
				contrl(e)
			}
		}
		contrl(nil)
		if sdl.GetTicks64()/1000 != secs {
			show = !show
		}
		secs = sdl.GetTicks64() / 1000
		renderer.SetDrawColor(0, 0, 0, 255)
		renderer.Clear()
		displist()
		loop()
		renderer.Present()
		sdl.Delay(1000 / 60)
	}

	renderer.Destroy()
	window.Destroy()
	sdl.Quit()
}
