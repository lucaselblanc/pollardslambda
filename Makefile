TARGET    := lambda
CXX       := g++

LDLIBS    := -pthread -ldl -lrt -lcrypto

SRC_CPP   := modinv.cpp lambda.cpp secp256k1.cpp
OBJ_CPP   := $(SRC_CPP:.cpp=.o)
OBJ       := $(OBJ_CPP)

.PHONY: all clean

all: $(TARGET)

CXXFLAGS  := -O3 -march=native -fomit-frame-pointer -std=c++14 -pthread -I. -MD -flto=auto -fno-rtti
LDFLAGS   := -fwhole-program -flto=auto

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	$(CXX) $(LDFLAGS) $(CXXFLAGS) $(OBJ) -o $@ $(LDLIBS)

clean:
	@echo "Cleaning..."
	rm -f $(TARGET)
	find . -type f \( -name "*.o" -o -name "*.d" \) -delete