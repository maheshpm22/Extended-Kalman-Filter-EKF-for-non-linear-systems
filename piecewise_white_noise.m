% Piece wise white noise function 

function output = piecewise_white_noise(ndim , variance , dt) 


if ndim == 5 

    G = @(t)[t^4/4 ; t^3/3 ; t^2/2 ; t ; 1 ] ; 
    Gc = @(t) [t^4/4  t^3/3  t^2/2 t  1 ] ; 
    output = variance*G(dt)*Gc(dt) ; 


elseif ndim == 4 

    G = @(t)[ t^3/3 ; t^2/2 ; t ; 1 ] ; 
    Gc = @(t) [ t^3/3  t^2/2 t  1 ] ; 
    output = variance*G(dt)*Gc(dt) ; 

elseif ndim ==3 

    G = @(t)[ t^2/2 ; t ; 1 ] ; 
    Gc = @(t) [  t^2/2 t  1 ] ; 
    output = variance*G(dt)*Gc(dt) ; 

elseif ndim == 2 

    G = @(t)[ t^2/2 ; t   ] ; 
    Gc = @(t) [  t^2/2 t  ] ; 
    output = variance*G(dt)*Gc(dt) ; 

end

end







