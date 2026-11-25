#set page(paper: "a4", margin: (x: 2cm, y: 1.5cm), numbering: "1")
#set text(lang: "pl", size: 10.5pt)
#set par(justify: true)
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)")
#show heading: set block(above: 1.4em, below: 0.8em)
// odwołanie do wzoru jako "(3)" zamiast "Równanie 3"
#show ref: it => {
  let el = it.element
  if el != none and el.func() == math.equation {
    link(el.location(), numbering(el.numbering, ..counter(math.equation).at(el.location())))
  } else {
    it
  }
}

#let wynik(body) = block(stroke: 1pt, inset: (x: 10pt, y: 8pt), radius: 3pt, width: 100%, body)
// (v Δt / Δx)^2 - tak jak na slajdach, bez skrótu
#let vdx = $((v Delta t)/(Delta x))^2$

#align(center)[
  #text(15pt, weight: "bold")[Równanie falowe z tłumieniem: warunki brzegowe odwrócone] \
  #v(2pt)
  Dirichlet na lewym końcu, Neumann na prawym końcu \
  (odpowiednik slajdów 33 i 35 z wykładu 5)
]

= Oznaczenia (jak na wykładzie)

- $y_m^((k)) = y(x_m, t_k)$, węzły $m = 1, dots, M$, poziomy czasu $k = 1, 2, dots$
- $beta = b/(2 mu)$, #h(1em) $g_m = g(x_m)$
- lewy koniec (Dirichlet, slajd 21): $y(x_1, t) = l(t)$, #h(0.5em) $l_k = l(t_k)$
- prawy koniec (Neumann, slajd 24): $(partial y)/(partial x) (x_M, t) = h(t)$, #h(0.5em) $h_k = h(t_k)$

= Wzór ogólny (slajd 29)

Równanie $(partial^2 y)/(partial t^2) = v^2 (partial^2 y)/(partial x^2) - 2 beta (partial y)/(partial t)$
po zastąpieniu pochodnych różnicami centralnymi:

$ (y_m^((k+1)) - 2 y_m^((k)) + y_m^((k-1)))/(Delta t^2)
  = v^2 (y_(m+1)^((k)) - 2 y_m^((k)) + y_(m-1)^((k)))/(Delta x^2)
  - beta (y_m^((k+1)) - y_m^((k-1)))/(Delta t) $

Mnożymy przez $Delta t^2$ i przenosimy $y_m^((k+1))$ na lewą stronę:

$ (1 + Delta t beta) y_m^((k+1))
  = vdx (y_(m+1)^((k)) - 2 y_m^((k)) + y_(m-1)^((k))) + 2 y_m^((k)) - (1 - Delta t beta) y_m^((k-1)) $ <ogolny>

Dla węzła brzegowego $m = M$ wzór potrzebuje sąsiada $y_(M+1)^((k))$, który leży *poza struną*
(węzeł-duch). Na wykładzie duch był po lewej ($y_0^((k))$), tutaj jest po prawej.

= Węzeł-duch z warunku Neumanna

Pochodną w $x_M$ zastępujemy różnicą centralną:

$ (y_(M+1)^((k)) - y_(M-1)^((k)))/(2 Delta x) = h_k
  quad => quad y_(M+1)^((k)) = y_(M-1)^((k)) + 2 Delta x thin h_k $ <duch>

Na slajdzie 25 jest $(y_2^((k)) - y_0^((k)))/(2 Delta x) = h_k$, skąd $y_0^((k)) = y_2^((k)) - 2 Delta x thin h_k$.
Tutaj znak jest *plus*, bo duch leży po drugiej stronie węzła brzegowego.

Wzór ogólny @ogolny zapisany dla węzła brzegowego $m = M$ (w nawiasie pojawia się duch):

$ (1 + Delta t beta) y_M^((k+1))
  = vdx (y_(M+1)^((k)) - 2 y_M^((k)) + y_(M-1)^((k))) + 2 y_M^((k)) - (1 - Delta t beta) y_M^((k-1)) $ <brzeg>

Na slajdzie 25 jest to samo dla lewego końca $m = 1$ i bez tłumienia ($beta = 0$):

#math.equation(block: true, numbering: none,
  $y_1^((k+1)) = vdx (y_2^((k)) - 2 y_1^((k)) + y_0^((k))) + 2 y_1^((k)) - y_1^((k-1))$)

Wstawiamy ducha @duch do nawiasu we wzorze @brzeg:

$ y_(M+1)^((k)) - 2 y_M^((k)) + y_(M-1)^((k))
  = 2 (y_(M-1)^((k)) - y_M^((k)) + Delta x thin h_k) $ <nawias>

= Prawy koniec, kolejne kroki ($k >= 2$)

Po wstawieniu @nawias do @brzeg i podzieleniu przez $(1 + Delta t beta)$:

#wynik[
$ y_M^((k+1)) = 2/(1 + Delta t beta) vdx (y_(M-1)^((k)) - y_M^((k)) + Delta x thin h_k)
  + 2/(1 + Delta t beta) y_M^((k))
  - (1 - Delta t beta)/(1 + Delta t beta) y_M^((k-1)) $ <krok-k>
]

= Prawy koniec, pierwszy krok ($k = 1$)

Wzór @krok-k dla $k = 1$ wymaga $y_M^((0))$ (chwila przed startem, drugi duch, tym razem w czasie).
Wyznaczamy go z prędkości początkowej:

$ (y_M^((2)) - y_M^((0)))/(2 Delta t) = g_M
  quad => quad y_M^((0)) = y_M^((2)) - 2 Delta t thin g_M $

Wstawiamy do @brzeg (dla $k = 1$) razem z @nawias:

$ (1 + Delta t beta) y_M^((2))
  = 2 vdx (y_(M-1)^((1)) - y_M^((1)) + Delta x thin h_1) + 2 y_M^((1))
  - (1 - Delta t beta)(y_M^((2)) - 2 Delta t thin g_M) $

Wyrazy z $y_M^((2))$ przenosimy na lewo, $(1 + Delta t beta) + (1 - Delta t beta) = 2$, i dzielimy przez 2:

#wynik[
$ y_M^((2)) = vdx (y_(M-1)^((1)) - y_M^((1)) + Delta x thin h_1)
  + y_M^((1)) + (1 - Delta t beta) Delta t thin g_M $ <krok-1>
]

= Lewy koniec (Dirichlet)

Wartość jest zadana wprost, nic nie trzeba wyprowadzać:

$ y_1^((k)) = l_k quad "dla każdego" k $

Węzły wewnętrzne $m = 2, dots, M-1$ liczy się bez zmian, tymi samymi wzorami co na wykładzie.

= Rozwiązanie równowagowe (slajd 35)

W stanie równowagi nic się nie zmienia w czasie:

$ (partial y)/(partial t) = (partial^2 y)/(partial t^2) = 0
  quad => quad (partial^2 y_E)/(partial x^2) = 0
  quad => quad y_E (x) = c x + d $

Stałe wyznaczamy z warunków brzegowych (teraz zamienionych miejscami):

$ y_E (0) = d = l, quad quad (partial y_E)/(partial x) (L) = c = h $

#wynik[
$ y_E (x) = h x + l $
]

Dla danych z zestawu 1 ($l = sqrt(2)\/2$, $h = -3 pi sqrt(2)\/2$):

$ y_E (x) = -(3 pi sqrt(2))/2 x + sqrt(2)/2 approx -6.6643 x + 0.7071 $

= Porównanie z wykładem

#table(
  columns: (auto, 1fr, 1fr),
  inset: 7pt,
  align: (left + horizon, center + horizon, center + horizon),
  table.header([], [*Wykład* \ Neumann lewy, Dirichlet prawy], [*Laboratorium* \ Dirichlet lewy, Neumann prawy]),
  [Warunki], [$(partial y)/(partial x)(x_1, t) = h(t), quad y(x_M, t) = r(t)$], [$y(x_1, t) = l(t), quad (partial y)/(partial x)(x_M, t) = h(t)$],
  [Węzeł-duch], [$y_0^((k)) = y_2^((k)) - 2 Delta x thin h_k$], [$y_(M+1)^((k)) = y_(M-1)^((k)) + 2 Delta x thin h_k$],
  [Nawias we wzorze], [$y_2^((k)) - y_1^((k)) - Delta x thin h_k$], [$y_(M-1)^((k)) - y_M^((k)) + Delta x thin h_k$],
  [Liczony węzeł], [$y_1$ (lewy)], [$y_M$ (prawy)],
  [Równowaga], [$c = h, quad d = r - h L$], [$c = h, quad d = l$],
)
