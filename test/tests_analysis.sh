
octave --eval runTestProblems_local > aux.txt

#more aux.txt | grep -wns "Solving problem:" -A 5 > aux.txt
