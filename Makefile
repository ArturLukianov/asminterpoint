SOURCES := main print memory strings input matrix method
OBJECTS := $(patsubst %, %.o, $(SOURCES))
FLAGS   := -no-pie -s -nostdlib --gc-sections -z noseparate-code
DFLAGS  := 

main: $(OBJECTS)
	ld -o main $(OBJECTS) $(DFLAGS)
%.o: %.S
	nasm -f elf64 $^