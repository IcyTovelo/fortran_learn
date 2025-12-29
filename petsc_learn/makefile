# 1. 基础路径配置
PETSC_DIR = /opt/homebrew/Cellar/petsc/3.24.2
EXTRA_LIB_PATH = -L/opt/homebrew/opt/openblas/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# ex1
ex1: ex1.o
	-${FLINKER} -o ex1 ex1.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex1.o

# ex2
ex2: ex2.o
	-${FLINKER} -o ex2 ex2.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex2.o

# ex3 (新增！矩阵测试)
ex3: ex3.o
	-${FLINKER} -o ex3 ex3.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex3.o

# ex4 (新增！KSP求解器)
ex4: ex4.o
	-${FLINKER} -o ex4 ex4.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex4.o

clean_files:
	${RM} ex1 ex2 ex3 ex4