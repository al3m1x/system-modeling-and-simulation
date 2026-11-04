#set page(paper: "a4", margin: 2.2cm, numbering: "1")
#set text(lang: "pl", size: 10pt)
#set par(justify: true)
#set math.equation(numbering: "(1)")
#show heading.where(level: 1): set text(size: 13pt)
#show heading.where(level: 2): set text(size: 11pt)
#show heading: set block(above: 1.4em, below: 0.9em)
#show table: set text(size: 9pt)

#let wyn(rows) = table(
  columns: (auto, auto, auto, auto, auto, auto, auto, auto),
  align: (center, right, right, right, right, right, right, right),
  stroke: none,
  inset: (x: 0.5em, y: 0.42em),
  table.hline(stroke: 0.5pt),
  table.header([\#], [$x_I$], [$y_I$], [$z_I$], [$Delta t$], [$E_"kin"$], [$E_"pot"$], [$E_"total"$]),
  table.hline(stroke: 0.3pt),
  ..rows.map(row => row.map(cell => [#cell])).flatten(),
  table.hline(stroke: 0.5pt),
)

#align(center)[
  #text(size: 11pt)[Modelowanie i symulacja systemów]

  #v(0.2em)
  #text(size: 14pt, weight: "bold")[Laboratorium 2: Odbijanie się ciała od powierzchni w jednorodnym polu grawitacyjnym]

  #v(0.6em)
  #text(size: 11pt)[
    Sebastian Kwaśniak, 188807 \
    Juliusz Radziszewski, 193504
  ]
]
#v(0.8em)

= Zadanie

Ciało o masie $m = 1 "kg"$ porusza się w jednorodnym polu grawitacyjnym ($g = 9.81 "m/s"^2$) i odbija się od powierzchni $z = f(x, y) = x^2 - y^2$. Warunki początkowe: $(x_0, y_0, z_0) = (-10, 0.01, 200)$ oraz $bold(v)_0 = (10, 0, 0)$. Symulujemy pięć pierwszych odbić, wyznaczając punkty odbicia, czasy $Delta t$ między kolejnymi kolizjami oraz energie $E_"kin"$, $E_"pot"$ i $E_"total"$ w punktach odbicia – dla odbić sprężystych ($k = 1$) oraz ze stratami energii ($k = 0.9$), gdzie $k$ oznacza stosunek energii kinetycznej po odbiciu do energii kinetycznej sprzed odbicia.

= Metoda rozwiązania

Ruch ciała opisują równania $x = x_0 + v_(0x) t$, $y = y_0 + v_(0y) t$ oraz $z = z_0 + v_(0z) t - 1/2 g t^2$. Po podstawieniu ich do warunku zetknięcia z powierzchnią $z(t) = f(x(t), y(t))$ otrzymujemy równanie kwadratowe $A t^2 + B t + C = 0$ (dla $f$ stopnia $<= 2$). Odrzucamy pierwiastki zespolone i niedodatnie; czas pierwszego odbicia $t_I$ to najmniejszy z pozostałych – przy pierwszym odbiciu drugi pierwiastek jest ujemny, a przy kolejnych bliski zeru, bo odpowiada punktowi startowemu.

Punkt odbicia wyznaczamy z równań ruchu dla $t = t_I$; prędkość tuż przed odbiciem wynosi $bold(v)^("before") = [v_(0x), v_(0y), v_(0z) - g t_I]$. Kierunek odbicia określa znormalizowany wektor normalny powierzchni w punkcie kolizji

$ bold(n) = bold(N)/norm(bold(N)), quad bold(N) = [-(partial f)/(partial x), -(partial f)/(partial y), 1], $

a prędkość po odbiciu wynosi

$ bold(v)^("after") = sqrt(k) [bold(v)^("before") - 2 (bold(v)^("before") dot bold(n)) bold(n)]. $

Punkt odbicia staje się nowym warunkiem początkowym i procedurę powtarzamy pięciokrotnie. Energie w punktach odbicia liczymy jako $E_"kin" = 1/2 m norm(bold(v))^2$ (z prędkością tuż przed odbiciem), $E_"pot" = m g z_I$ oraz $E_"total" = E_"kin" + E_"pot"$. Program (Python, NumPy + Matplotlib) rozwiązuje równanie kolizji funkcją `np.roots`, a trajektorie rysuje na siatce powierzchni.

= Wyniki

W tabelach pierwszy wiersz podaje czas od startu do pierwszego odbicia, kolejne – czasy między sąsiednimi kolizjami; $Delta t$ w sekundach, energie w dżulach.

== Odbicia sprężyste ($k = 1$)

#figure(
  align(center, wyn((
    ("1", "13.1776", "0.0100", "173.6502", "2.3178", "308.4917", "1703.5083", "2012.0000"),
    ("2", "−10.4561", "0.0433", "109.3283", "2.0210", "939.4889", "1072.5111", "2012.0000"),
    ("3", "7.5090", "0.1922", "56.3475", "1.1499", "1459.2308", "552.7692", "2012.0000"),
    ("4", "−5.2287", "0.8205", "26.6657", "0.5709", "1750.4096", "261.5904", "2012.0000"),
    ("5", "5.0232", "3.9680", "9.4880", "0.3344", "1918.9231", "93.0769", "2012.0000"),
  ))),
  kind: table,
  caption: [Punkty odbicia, czasy między kolizjami i energie dla odbić sprężystych ($k = 1$).],
)

Całkowita energia jest stała: $E_"total" = 2012.0000 "J" = m g z_0 + 1/2 m norm(bold(v)_0)^2$. Energia kinetyczna rośnie, a potencjalna maleje – zamieniają się one wyłącznie między sobą. Zwraca uwagę naprzemienny znak $x_I$: ciało odbija się na przemian po obu stronach siodła, a czas lotu skraca się z każdym odbiciem (z 2.32 s do 0.33 s).

#figure(image("fig_bounces_k1.0.png", width: 82%), kind: image, caption: [Pięć pierwszych odbić od powierzchni $z = x^2 - y^2$ bez strat energii ($k = 1$); punkty odbić zaznaczono czarnymi kropkami.])

== Straty energii ($k = 0.9$)

#figure(
  align(center, wyn((
    ("1", "13.1776", "0.0100", "173.6502", "2.3178", "308.4917", "1703.5083", "2012.0000"),
    ("2", "−10.3655", "0.0431", "107.4424", "2.1221", "927.1409", "1054.0099", "1981.1508"),
    ("3", "7.2617", "0.1889", "52.6966", "1.2349", "1371.4836", "516.9532", "1888.4367"),
    ("4", "−4.8780", "0.7874", "23.1748", "0.6089", "1523.9431", "227.3453", "1751.2884"),
    ("5", "4.5674", "3.7090", "7.1042", "0.3505", "1529.2016", "69.6924", "1598.8941"),
  ))),
  kind: table,
  caption: [Punkty odbicia, czasy między kolizjami i energie dla $k = 0.9$.],
)

Pierwsze odbicie jest identyczne jak w przypadku sprężystym. Od drugiego odbicia całkowita energia maleje, ponieważ po odbiciu zachowana zostaje energia potencjalna $E_"pot" (z_I)$ oraz część $k$ energii kinetycznej sprzed odbicia (np. $1703.5083 + 0.9 dot 308.4917 = 1981.1508$). Kolejne wartości $E_"total"$ maleją w przybliżeniu geometrycznie, a czasy między kolizjami skracają się.

#figure(image("fig_bounces_k0.9.png", width: 82%), kind: image, caption: [Pięć pierwszych odbić od powierzchni $z = x^2 - y^2$ ze stratami energii ($k = 0.9$).])

= Wnioski

+ Warunek kolizji $z(t) = f(x(t), y(t))$ sprowadza się do równania kwadratowego względem $t$; najmniejszy dodatni pierwiastek rzeczywisty daje czas pierwszego odbicia, również dla kolejnych odbić.
+ Przy odbiciach sprężystych ($k = 1$) całkowita energia mechaniczna jest zachowana (2012 J w każdym punkcie odbicia); zmienia się jedynie podział na energię kinetyczną i potencjalną.
+ Dla $k = 0.9$ całkowita energia maleje z każdym odbiciem (z 2012 J do 1598.89 J po pięciu odbiciach), a czasy między kolizjami skracają się (z 2.32 s do 0.35 s) – ciało stopniowo traci energię.
