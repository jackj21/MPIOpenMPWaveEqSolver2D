INCLUDE_DIR := include
SRC_DIR := src

OBJECTS := array_2d.o wave.o

C_FLAGS := 

INCLUDES := -lm

all: wave_timing wave_images wave_animation wave_error

array_2d.o:
	mpicc -o array_2d.o -I$(INCLUDE_DIR) -c $(SRC_DIR)/array_2d.c $(C_FLAGS) $(INCLUDES)

wave.o: array_2d.o
	mpicc -o wave.o -I$(INCLUDE_DIR) -c $(SRC_DIR)/wave.c $(C_FLAGS) $(INCLUDES)

wave_timing: $(OBJECTS)
	mpicc -o wave_timing wave_timing.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_images: $(OBJECTS)
	mpicc -o wave_images wave_images.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_animation: $(OBJECTS)
	mpicc -o wave_animation wave_animation.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

wave_error: $(OBJECTS)
	mpicc -o wave_error wave_error.c -I$(INCLUDE_DIR) $(OBJECTS) $(C_FLAGS) $(INCLUDES)

clean:
	rm -f *.o
	rm -f wave_timing
	rm -f wave_images
	rm -f wave_animation
	rm -f wave_error
