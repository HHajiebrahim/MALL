function [f,BestAnswer,Error,LR,UR,Kind,Name,Step,nVar]=f_3(x)
    
    f=x(1)^2+1e6*sumsqr(x(2):x(end));

    BestAnswer = 0;
    
    Error = abs (f-BestAnswer);
    
    LR= -100;
    UR= 100;
    
    Kind='min';
    
    Name='Cigar';
    
    Step=2.5;

    nVar = 30;

end