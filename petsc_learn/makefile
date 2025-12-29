# 1. 基础路径配置
PETSC_DIR = /opt/homebrew/Cellar/petsc/3.24.2
EXTRA_LIB_PATH = -L/opt/homebrew/opt/openblas/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# 2. 编译规则
# ex1 (之前的)
ex1: ex1.o
	-${FLINKER} -o ex1 ex1.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex1.o

# ex2 (新增的！)
ex2: ex2.o
	-${FLINKER} -o ex2 ex2.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex2.o

# 清理命令
clean_files:
	${RM} ex1 ex2