# ======================(A)
Only in the cmd (choose the Command Prompt )
(1) "R" and "R --version" can start
(2) "postitron ." can start postitron program

in Powersehll,
(1) "R" or "r" 是 PowerShell 的内置别名，用于调用命令历史.  AND "R" or "R --version" has no reponse
(2) "positron ." produce show error message

# ======================(B)
Current activate within Conda 

Positon App works with different Python very well
Positon cannnot detect conda R at all, only the system R 4.5.1 are activated everytime



C:\Users\zhen-\anaconda3\Lib\R\bin>conda env list

# conda environments:
#
                       C:\Users\zhen-\Anaconda3
DSp3.10                C:\Users\zhen-\Anaconda3\envs\DSp3.10
R43                    C:\Users\zhen-\Anaconda3\envs\R43
d2l                    C:\Users\zhen-\Anaconda3\envs\d2l
sc_omics               C:\Users\zhen-\Anaconda3\envs\sc_omics
torch_cuda             C:\Users\zhen-\Anaconda3\envs\torch_cuda
base                   C:\Users\zhen-\anaconda3
DSp3.10                C:\Users\zhen-\anaconda3\envs\DSp3.10
R43                    C:\Users\zhen-\anaconda3\envs\R43
d2l                    C:\Users\zhen-\anaconda3\envs\d2l
sc_omics               C:\Users\zhen-\anaconda3\envs\sc_omics
torch_cuda             C:\Users\zhen-\anaconda3\envs\torch_cuda
DSp3.10                c:\Users\zhen-\anaconda3\envs\DSp3.10
d2l                    c:\Users\zhen-\anaconda3\envs\d2l
sc_omics               c:\Users\zhen-\anaconda3\envs\sc_omics
torch_cuda             c:\Users\zhen-\anaconda3\envs\torch_cuda


(R43) C:\Users\zhen->where R

C:\Users\zhen-\anaconda3\envs\R43\Scripts\R.exe
C:\Users\zhen-\Anaconda3\Scripts\R.exe
C:\Program Files\R\R-4.5.1\bin\R.exe


# positoron default terminal is powershell

1) cmd  #changet to cmd 
2) conda activate R43
3) R    # ====> the topright can see the R 4.3. 3 very Clearly!
4) Then Ctrl + Enter to run each line in file
5) q() # or quit() retrun to cmd 

1A) cmd  #changet to cmd 
2A) conda activate R41
3A) R # ====> the topright can see the R 4.2. 3 very Clearly!
4A) Then Ctrl + Enter to run each line in file
5A) q(A) # or quit(A) retrun to cmd 

1) cmd  #changet to cmd 
2) conda activate R43
3) R.exe    # something R it self is a alias of other things 
4) Then Ctrl + Enter to run each line in file
5) q() # or quit() retrun to cmd 

1A) cmd  #changet to cmd 
2A) conda activate R41
3A) R # ====> the topright can see the R 4.2. 3 very Clearly!
4A) Then Ctrl + Enter to run each line in file
5A) q(A) # or quit(A) retrun to cmd 

