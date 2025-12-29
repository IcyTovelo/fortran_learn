
# 1. PETSc 路径
PETSC_DIR = /opt/homebrew/Cellar/petsc/3.24.2

# 2. 告诉链接器库文件在哪里
# (A) /opt/homebrew/opt/openblas/lib  -> 解决 lapack not found
# (B) /Library/.../MacOSX.sdk/usr/lib -> 解决 c++ not found
EXTRA_LIB_PATH = -L/opt/homebrew/opt/openblas/lib -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

# 3. 引入 PETSc 标准规则
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

# 4. 编译目标
ex1: ex1.o
	-${FLINKER} -o ex1 ex1.o ${PETSC_LIB} ${EXTRA_LIB_PATH}
	${RM} ex1.o

clean_files:
	${RM} ex1