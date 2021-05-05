INCLUDE_DIR := include
SRC_DIR := src

OBJECTS := array_2d.o wave.o array_2d_io.o

C_FLAGS := 

INCLUDES := -lm -fopenmp

obj: array_2d.o wave.o array_2d_io.o

all: wave_timing wave_images wave_animation wave_error wave_print_test

array_2d.o:
	mpicc -o array_2d.o -I$(INCLUDE_DIR) -c $(SRC_DIR)/array_2d.c $(C_FLAGS) $(INCLUDES)

wave.o: array_2d.o
	mpicc -o wave.o -I$(INCLUDE_DIR) -c $(SRC_DIR)/wave.c $(C_FLAGS) $(INCLUDES)

array_2d_io.o: 
	mpicc -o array_2d_io.o -I$(INCLUDE_DIR) -c $(SRC_DIR)/array_2d_io.c $(C_FLAGS) $(INCLUDES)

wave_timing: $(OBJECTS)
	mpicc -o wave_timing wave_timing.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_images: $(OBJECTS)
	mpicc -o wave_images wave_images.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_animation: $(OBJECTS)
	mpicc -o wave_animation wave_animation.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_error: $(OBJECTS)
	mpicc -o wave_error wave_error.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_print_test: $(OBJECTS)
	mpicc -o wave_print_test wave_print_test.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)



clean:
	rm -f *.o
	rm -f wave_timing
	rm -f wave_images
	rm -f wave_animation
	rm -f wave_error
	rm -f wave_print_test
