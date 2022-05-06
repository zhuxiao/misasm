CXXFLAGS = -O0 -g -Wall -fmessage-length=0 -std=c++11

OBJS = alnDataLoader.o Base.o blatAlnTra.o Block.o Chrome.o clipAlnDataLoader.o clipReg.o covLoader.o events.o FastaSeqLoader.o Genome.o InsertSizeEst.o LocalAssembly.o misAlnReg.o misasm_main.o Paras.o RefSeqLoader.o Region.o Thread.o util_header.o util.o varCand.o 

LIBS += -lhts -lpthread   

TARGET = misasm

all: $(TARGET) clean

$(TARGET): $(OBJS)
	$(CXX) -o $(TARGET) $(OBJS) $(LIBS) 

clean:
	rm -f $(OBJS)
	
clean-all: clean
	rm -f $(TARGET)
