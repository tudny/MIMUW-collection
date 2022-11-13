
# Moduł ROBSON
Aby korzystać z modułu ROBSON, należy skompilować go poniższą komendą.
```bash
mvn clean compile
```
Aby wygenerować plik ```.jar``` wystarczy komenda
```bash
mvn clean compile assembly:single
```
Wówczas zostanie stworzony plik zawierający cały moduł wraz z potrzebnymi bibliotekami.
Znajduje się on w katalogu ```target``` i jego nazwa to ```ROBSON-wersja-jar-with-dependencies.jar```.
## Dostępne komendy
Dostępne są wszystkie komendy zadane według treści zadania.
## fromJson()
Wedle treści zadania
## toJson()
Wedle treści zadania
## toJava()
Wedle treści zadania.
Wygenerowany program można skompilować i uruchomić za pomocą
```bash
javac nazwa.java
java nazwa
```
## Testowanie
Moduł ROBSON wyposażony jest w pewien zestaw testów. Aby uruchomić wszystkie testy jednostkowe
należy wykonać polecenie
```bash
mvn test
```

# Moduł Henry
## Uruchomienie
Henry to przyjazny ROB żyjący na małej planszy potrzebujący instrukcji i baterii, aby przeżyć.
Aby rozpocząć przygodę z Henrym, należy
- uruchomić program za pomocą komendy
```bash
mvn clean javafx:run
```
- lub skompilować moduł ROBSON i uruchomić powstały plik z rozszerzeniem ```.jar``` za pomocą
```bash
java -jar ROBSON-wersja-jar-with-dependencies.jar
```
### Parametry planszy
Program testowy był na monitorze ```1920x1080```, zatem nie musi być zawsze dobrze wyświetlany.
Domyślne ustawienia zapisane są w pliku konfiguracyjnym ```src/main/java/pl.edu.mimuw.at429630.henry.BoardConfig```.
Domyślne wartości można nadpisać na czas uruchomienia programu, dodając argumenty podczas uruchamiania.
Przykładowe użycie (nie trzeba uzupełnić wszystkich parametrów, kolejność dowolna).
```
  java -jar ROBSON-wersja-jar-with-dependencies.jar -TW100 -TW100 -BW8 -BH6
```
Wyjaśnienie:
- -TW100 (tail width 100) ustawia szerokość pola na 100 pikseli.
- -TH100 (tail height 100) ustawia wysokość pola na 100 pikseli.
- -BW8 (board width 8) ustawia szerokość planszy na 8 pól.
- -BH6 (board height 6) ustawia wysokość planszy na 6 pól.


## Użycie
Klikając na poszczególne pola na planszy, możemy dodawać baterie (pożywienie) potrzebne Henry'emu.
Na dole ekranu znajdują się cztery istotne elementy interfejsu użytkownika. Idąc od lewej jest to
- pole z informacjami
- stan energii Henry'ego
- pole poleceń
- zatwierdzenie poleceń (Można zatwierdzać poprzez ENTER.)

Dostępne komendy to
- ```LOAD <ścieżka do pliku>``` - wczytanie programu w języku ROBSON
- ```RUN``` - uruchomienie programu 

Rozszerzenie języka ROBSON pomocne w sterowaniu Henrym
- ```Lewo```
- ```Prawo```
- ```Prosto```
- ```Jedz```
- ```Wachaj```

Należy używać ich na przykład w następujący sposób w pliku ```JSON```.
```json
{
  "typ": "Lewo"
}
```
## Przykładowy program sterujący Henrym
```json
{
  "typ": "Blok",
  "instrukcje": [
    {
      "typ": "Lewo"
    },
    {
      "typ": "Prosto"
    },
    {
      "typ": "Prawo"
    },
    {
      "typ": "Jedz"
    },
    {
      "typ": "Wachaj"
    }
  ]
}

```

## Przykładowy scenariusz
- Uruchamiamy program okienkowy.
- Ustawiamy kilka baterii na planszy poprzez klikanie na pola.
- Wpisujemy ```LOAD examples/henry.json```.
- Uruchamiamy program przez ```RUN```.
- (Ten sam program można uruchamiać wiele razy, bez ponownego wczytywania.)
- Henry wczytuje cały program i wykonuje go instrukcje po instrukcji co sekundę. Cały program uruchamiany jest od razu,
  lecz polecenia ruchu wykonywane są później, aby mogły być widoczne dla użytkownika.
- Jeżeli Henry'emu skończy się energia, umrze.

Film z podanym scenariuszem https://youtu.be/8QiAZSurh6Y

# Źródła
- Grafika trawy: (`src\main\resources\textures\terrain_texture.jpg`) https://freestocktextures.com/texture/green-grass-background-2,944.html
- Grafika z ramką: powyższa grafika z autorską modyfikacją
- Grafika Henry'ego: autorska grafika
- Grafika Baterii: autorska grafika

# Błędy i zastrzeżenia
- Program testowany na systemie Windows 10.
- Używana wersja JRE, JDK, JavaFX:
```bash
java version "11.0.11" 2021-04-20 LTS
Java(TM) SE Runtime Environment 18.9 (build 11.0.11+9-LTS-194)
Java HotSpot(TM) 64-Bit Server VM 18.9 (build 11.0.11+9-LTS-194, mixed mode)
```
```bash
javac 11.0.11
```
```bash
openjfx 11
```
