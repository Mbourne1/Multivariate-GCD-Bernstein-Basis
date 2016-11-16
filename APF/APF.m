function [fxy_lr, gxy_lr, uxy_lr, vxy_lr, dxy_lr, alpha_lr, th1_lr, th2_lr] = ...
    APF(fxy,gxy,uxy,vxy,m,n,t,t1,t2,alpha,th1,th2)

global SETTINGS
switch SETTINGS.APF_METHOD
    
    
    case 'Standard APF Nonlinear'
        
        error('err');
        
    case 'None'
        
        fww = GetWithThetas(fxy,th1,th2);
        a_gww = alpha.*GetWithThetas(gxy,th1,th2);
        uww = GetWithThetas(uxy,th1,th2);
        vww = GetWithThetas(vxy,th1,th2);
        
        [dww] = GetGCD_Coefficients(fww,a_gww,uww,vww,m,n,t,t1,t2);
        
        dxy_lr = GetWithoutThetas(dww,th1,th2);
        
        fxy_lr = fxy;
        gxy_lr = gxy;
        uxy_lr = uxy;
        vxy_lr = vxy;
        alpha_lr = alpha;
        th1_lr = th1;
        th2_lr = th2;
        
end


end