#set document(
  title: "Laboratorium 1: Ruch ciała w jednorodnym polu grawitacyjnym",
  author: ("Sebastian Kwaśniak", "Juliusz Radziszewski"),
)
#set text(lang: "pl", size: 11pt)
#set par(justify: true, leading: 0.65em)
#set page(margin: 2.2cm, numbering: "1", number-align: center)
#set figure(gap: 0.8em)

#align(center)[
  #text(size: 16pt, weight: "bold")[Modelowanie i symulacja systemów]
  #v(0.3em)
  #text(size: 13pt)[Laboratorium 1: Ruch ciała w jednorodnym polu grawitacyjnym]
  #v(1em)
  Sebastian Kwaśniak, 188807 \
  Juliusz Radziszewski, 193504
]

= Zadanie

Ciało startuje z punktu $(0, 0)$ z prędkością
początkową $v_0$ skierowaną pod kątem $alpha$. Zadanie polega na znalezieniu takich par $(v_0, alpha)$ z obszaru
$v_0 in [1, 30] "m/s"$, $alpha in [1 degree, 89 degree]$, aby trajektoria
przecinała oba regiony dozwolone - $R_1$ i $R_2$ (odcinki poziome na wysokości
$y = 3$: $x in [1, 5]$ oraz $x in [9, 13]$) i nie przecinała
regionów zabronionych $R_3$ ani $R_4$ (odcinki pionowe $x = 7$: $y in [4, 5]$
oraz $y in [7, 8]$). Pary $(v_0, alpha)$ losowano
$N = "100 000"$ razy.

= Metoda rozwiązania

Zamiast całkowania ruchu z krokiem czasowym zastosowano analityczny test
przecięcia złożony z dwóch kroków: wyznaczenia chwil przecięcia oraz
sprawdzenia tych chwil.

*Wyznaczenie chwil przecięcia* (funkcja `crosses_segment`). Punkty odcinka
opisuje parametr $s in [0, 1]$: $P(s) = (x_1 + s d_x, y_1 + s d_y)$, gdzie
$d_x = x_2 - x_1$, $d_y = y_2 - y_1$. Podstawiając równania ruchu do prostej
zawierającej odcinek i eliminując $s$, otrzymujemy równanie kwadratowe
względem czasu:

$ A t^2 + B t + C = 0 quad "gdzie" quad A = -1/2 d_x g, quad B = d_x v_(0y) - d_y v_(0x), quad C = d_y x_1 - d_x y_1 $

Dla odcinków pionowych ($d_x = 0$), takich jak $R_3$ i $R_4$, równanie jest
liniowe ($t = -C \/ B$). Jeżeli $Delta = B^2 - 4 A C >= 0$, rozwiązania
$t_1$ i $t_2$ to chwile, w których trajektoria spotyka prostą zawierającą
odcinek.

*Test chwil* (funkcja `hit_at_time`). Każdą z tych chwil sprawdzamy
niezależnie: przecięcie musi nastąpić w trakcie lotu ($0 <= t <= t_f$, gdzie
$t_f = 2 v_(0y) \/ g$ to czas powrotu na poziom $y = 0$), a ciało musi leżeć
na odcinku, a nie jedynie na jego przedłużeniu - w implementacji warunek ten
wyraża parametr $s$ wyznaczony z położenia ciała, który musi spełniać
$s in [0, 1]$. Odcinek jest przecięty, gdy którykolwiek z pierwiastków
spełnia oba warunki, a region - gdy przecięty jest dowolny jego odcinek.
Wszystkie sprawdzenia wykonano wektorowo w NumPy.

= Wyniki

W jednym z uruchomień programu (wygenerowało ono również oba poniższe wykresy)
warunki spełniło $n = 1345$ par, co daje oszacowanie
$hat(p) = n \/ N approx 0.0135$.

#figure(
  image("fig_initial_values.png", width: 70%),
  caption: [
    Pary $(alpha, v_0)$ spełniające warunki (czarne punkty) na tle
    przeszukiwanego obszaru (ramka). Zbiór tworzy trzy pasma, odpowiadające
    przelotowi pod $R_3$ ($alpha approx [40 degree, 54 degree]$), przez
    szczelinę między $R_3$ i $R_4$ ($alpha approx [52 degree, 72 degree]$)
    oraz nad $R_4$ ($alpha approx [66 degree, 73 degree]$). We wszystkich
    pasmach $v_0 in [11, 16] "m/s"$.
  ],
)

#figure(
  image("fig_trajectories.png", width: 70%),
  caption: [
    Regiony dozwolone $R_1$, $R_2$ (zielone) i zabronione $R_3$, $R_4$
    (czerwone) oraz pięć trajektorii spełniających warunki (parametry
    początkowe w legendzie). Każda startuje z $(0, 0)$, przecina $R_1$ i $R_2$
    oraz omija przeszkody w przekroju $x = 7$.
  ],
)

= Dokładność oszacowania

Oszacowana częstość $hat(p) approx 1.35%$ ma błąd standardowy
$sigma = sqrt(hat(p)(1 - hat(p)) \/ N) approx 3.6 dot 10^(-4)$, czyli ok. 2.7%
wartości $hat(p)$. Dokładność zależy przede wszystkim od liczby prób $N$ -
błąd maleje jak $1 \/ sqrt(N)$ - oraz od wartości $p$ (przez czynnik
$p(1 - p)$). Przy jednorodnym i niezależnym losowaniu jedynym źródłem
niepewności jest skończona liczba prób.
