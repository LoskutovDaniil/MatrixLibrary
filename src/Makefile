FLAGS = -Wall -Wextra -Werror -std=c11
LFLAGS = -lcheck -lpthread -lm
SOURCES=s21_*.c 
OBJECTS=s21_*.o

ifeq ($(shell uname -s), Linux)
	LFLAGS = -lpthread -lcheck -lsubunit -lm 
endif

all: clean s21_matrix.a

# для сборки всех .c файлов в исполняемы
test: clean s21_matrix.a
	gcc $(FLAGS) -c test.c 
	gcc $(FLAGS) test.o $(OBJECTS) $(LFLAGS) -L. s21_matrix.a -o test.out
	rm -rf *.o 
	clear
	./test.out 

# для отладки
rebuild: clean
#	gcc $(FLAGS) -c
	gcc $(FLAGS) $(SOURCES) test.c  $(LFLAGS) -o a.out
	rm -rf *.o 

# Сборка библиотеки
s21_matrix.a:
	gcc $(FLAGS) -c $(SOURCES)
	ar rc s21_matrix.a $(OBJECTS)
	ranlib s21_matrix.a
#	ranlib не обязателен, не понятно кто и зачем его начал использовать

gcov_report: clean add_coverage test
	./test.out
	lcov -t "test" -o test.info -c -d .
	genhtml -o report test.info
	open report/index.html

checks:
	clang-format -n *.c *.h
	clang-format -style="{BasedOnStyle: Google}" -i *.c
	clang-format -style="{BasedOnStyle: Google}" -i *.h
	cppcheck enable=all missingIncludeSystem *.c

clean:
	rm -rf *.o *.out gcov_report *.a *.gcda *.gcno test test.info s21_matrix.a

add:
	git add Makefile s21_*.c s21_matrix.h .clang-format

.PHONY: commit

commit:
	@echo "Введите сообщение для коммита:"
	@read commit_message; \
	git add Makefile s21_*.c s21_matrix.h .clang-format test.c; \
	git commit -m "$$commit_message"; \
	git push origin develop

add_coverage:
	$(eval FLAGS += --coverage)