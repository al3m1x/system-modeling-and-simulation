#set page(paper: "a4", margin: 2.2cm, numbering: "1")
#set text(lang: "pl", size: 10pt)
#set par(justify: true)
#set math.equation(numbering: "(1)", supplement: none)
#show heading.where(level: 1): set text(size: 13pt)
#show heading.where(level: 2): set text(size: 11pt)
#show heading: set block(above: 1.4em, below: 0.9em)
#show table: set text(size: 9pt)

#align(center)[
  #text(size: 11pt)[Modelowanie i symulacja systemów]

  #v(0.2em)
  #text(size: 14pt, weight: "bold")[Laboratorium 3: Równanie falowe - drgania struny]

  #v(0.6em)
  #text(size: 11pt)[
    Sebastian Kwaśniak, 188807 \
    Juliusz Radziszewski, 193504
  ]
]
#v(0.8em)

= Zadanie

Symulujemy drgania struny opisane równaniem falowym z tłumieniem

$ (partial^2 y)/(partial t^2) = v^2 (partial^2 y)/(partial x^2) - 2 beta (partial y)/(partial t), $ <falowe>

gdzie $v = sqrt(T \/ mu)$ to prędkość fali, a $beta$ - współczynnik wygasania drgań. Dane: gęstość liniowa $mu = 0.01 "kg/m"$, naciąg $T = 100 "N"$, długość $L = 1 "m"$, $M = 100$ węzłów w przestrzeni i $N = 2000$ poziomów czasowych. Warunki początkowe i brzegowe:

$ y(x, 0) = sin(3 pi x + pi/4), quad (partial y)/(partial t)(x, 0) = x - x^2, quad y(0, t) = sqrt(2)/2, quad (partial y)/(partial x)(L, t) = -(3 pi sqrt(2))/2 approx -6.664. $

Lewy koniec struny ma zadane wychylenie (warunek Dirichleta), a prawy - nachylenie (warunek Neumanna). Warunki brzegowe są zgodne z kształtem początkowym: $s(0) = sin(pi\/4) = sqrt(2)\/2$ oraz $s'(L) = 3 pi cos(13 pi \/ 4) = -3 pi sqrt(2) \/ 2$. Zadanie obejmuje symulację bez tłumienia ($b = 0$), z tłumieniem ($b = 0.5$) oraz wyznaczenie rozwiązania równowagowego i porównanie go z ostatnią iteracją symulacji z tłumieniem.

= Metoda rozwiązania

Zastosowano jawny schemat różnic skończonych. Dziedzinę dzielimy na regularną siatkę $x_m = m Delta x$, $t_k = k Delta t$ i oznaczamy $y_m^k = y(x_m, t_k)$. Ponieważ $M$ to liczba węzłów, $Delta x = L \/ (M - 1) approx 0.0101 "m"$. Prędkość fali wynosi $v = sqrt(100 \/ 0.01) = 100 "m/s"$, a krok czasowy przyjęto zgodnie z poleceniem jako $Delta t = Delta x \/ v approx 1.01 dot 10^(-4) "s"$ - jest to największy krok spełniający warunek stabilności CFL ($v Delta t \/ Delta x <= 1$). Cała symulacja obejmuje zatem $t in [0, 0.202] "s"$.

Współczynnik tłumienia $b$ (tarcie na jednostkę długości) przeliczamy na $beta$ z porównania $b \/ mu = 2 beta$, czyli $beta = b \/ (2 mu)$; dla $b = 0.5$ daje to $beta = 25 "s"^(-1)$.

Wyprowadzenie schematu różnicowego - wzór ogólny, pierwszy krok czasowy, warunek Dirichleta na lewym końcu i warunek Neumanna (węzeł fikcyjny) na prawym końcu - oraz rozwiązanie równowagowe przedstawiono na rysunku 1.

#figure(image("rysunek1.jpg", width: 85%), kind: image, caption: [Wyprowadzenie schematu różnicowego (cz. 1).]) <kartka1>

Program napisano w Pythonie (NumPy + Matplotlib); w każdym kroku czasowym wszystkie węzły wewnętrzne obliczane są wektorowo. Wyniki $y(x, t)$ przedstawiono jako powierzchnie 3D (`plot_surface`).

= Wyniki

== Drgania bez tłumienia ($b = 0$)

#figure(image("fig_surface_b0.0.png", width: 78%), kind: image, caption: [Rozwiązanie $y(x, t)$ bez tłumienia ($b = 0$).])

Struna drga okresowo z niezmienną amplitudą. Prawy (swobodny) koniec przemieszcza się między $-0.40 "m"$ a $-11.51 "m"$, czyli symetrycznie wokół $-5.96 "m"$ (wartości równowagowej w $x = L$), a lewy koniec pozostaje nieruchomy w $y = sqrt(2)\/2 approx 0.707 "m"$. Kolejne maksima wychylenia prawego końca występują co $0.04 "s"$ ($t approx 0.0016, 0.0416, 0.0816, ... "s"$) i mają tę samą wartość. W czasie symulacji ($0.202 "s"$) widać zatem pięć pełnych okresów drgań.

== Drgania z tłumieniem ($b = 0.5$)

#figure(image("fig_surface_b0.5.png", width: 78%), kind: image, caption: [Rozwiązanie $y(x, t)$ z tłumieniem ($b = 0.5$, $beta = 25 "s"^(-1)$).])

Z tłumieniem amplituda drgań maleje z czasem. Najniższe wychylenie prawego końca wynosi już tylko $-9.34 "m"$ (wobec $-11.51 "m"$ bez tłumienia). Odchylenie prawego końca od równowagi w kolejnych maksimach wynosi $5.55$, $2.06$, $0.76$, $0.28$ i $0.10 "m"$ - każde jest ok. $0.37$ razy mniejsze od poprzedniego. Po kilku okresach struna praktycznie zastyga w położeniu prostoliniowym.

== Rozwiązanie równowagowe

Wyznaczenie rozwiązania równowagowego przedstawiono na rysunku 1.

Otrzymane rozwiązanie: $y_E (x) approx -6.664 x + 0.707$.

#figure(image("fig_equilibrium_b0.5.png", width: 66%), kind: image, caption: [Rozwiązanie równowagowe $y_E (x)$ (czerwona linia) i ostatnia iteracja symulacji z tłumieniem $b = 0.5$ (niebieska linia).])

Ostatnia iteracja symulacji niemal pokrywa się z rozwiązaniem równowagowym - maksymalna różnica $max_x |y(x, t_"end") - y_E (x)|$ wynosi $0.0371 "m"$ i występuje na prawym (swobodnym) końcu. Na początku symulacji odchylenie od równowagi $s(x) - y_E (x)$ wynosiło maksymalnie ok. $5.55 "m"$, więc do końca symulacji zmalało ok. 150-krotnie. Dla porównania bez tłumienia różnica ta w ostatniej iteracji wynosi $5.53 "m"$ - struna nadal drga z pełną amplitudą.

= Wnioski

+ Jawny schemat różnic skończonych z krokiem $Delta t = Delta x \/ v$ (granica warunku CFL) jest stabilny i dla $b = 0$ struna drga okresowo (okres $0.04 "s"$) ze stałą amplitudą.
+ Warunek brzegowy Neumanna realizujemy przez węzeł fikcyjny, a warunek prędkości początkowej - przez eliminację wartości $y^(-1)$ w pierwszym kroku; oba przypadki sprowadzają się do zmodyfikowanej postaci wzoru ogólnego.
+ Tłumienie ($beta = 25 "s"^(-1)$) powoduje zanik drgań i zbieżność rozwiązania do stanu równowagi $y_E (x) approx -6.664 x + 0.707$; ostatnia iteracja różni się od niego o najwyżej $0.0371 "m"$.
