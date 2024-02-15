PROLOG = swipl
INDEKS = at429630
ENTRYPOINT = program
BINARY_NAME = $(ENTRYPOINT)

all: $(BINARY_NAME)

$(BINARY_NAME): $(INDEKS).pl
	@echo "Compiling $<..."
	$(PROLOG) --goal=$(ENTRYPOINT) --stand_alone=true -o $@ -c $<

run: $(BINARY_NAME) tatry.txt
	@echo "Running $<..."
	@./$(BINARY_NAME) tatry.txt

clean:
	rm -f $(BINARY_NAME)

verify: verify.py $(INDEKS).pl
	./verify.py $(INDEKS).pl

convert: to_students.py $(INDEKS).pl
	./to_students.py $(INDEKS).pl
