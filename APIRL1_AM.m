

function out =  APIRL1_AM(f, Img, K, opts, optsNcvxAcc)



rho  = opts.rho;
tol = opts.tol;
Nit = opts.Nit;
p    = opts.p;

mu   = optsNcvxAcc.mu;
beta = optsNcvxAcc.beta;

err = 0;


relError        = zeros(Nit,1);
psnrGain        = relError;     % PSNR improvement every iteration
ssimGain        = relError;
objVal          = relError;

[row, col]  = size(f);

x = f;

%%%%%%% Acceleration Parameters %%%%%%%%
t = 1;
f = x;
y = x;

%*************** v sub-problem variable initialization ******************

v_1         = zeros(row, col); % v solution of the v sub-problem for Dx_h
v_2         = v_1; % v solutiono for the v sub-problem for Dx_v
v           = v_1;

% A is the blur kernel K, A = K;
eigA        = psf2otf(K,[row col]); %In the fourier domain
eigAtA      = abs(eigA).^2;
eigDtD      = abs(fft2([1 -1], row, col)).^2 + abs(fft2([1 -1]', row, col)).^2;


[D,Dt]      = defDDt(); %Declare forward finite difference operators

[Dx1, Dx2] = D(x);

lhs     = beta*eigDtD + eigAtA ;

tg = tic;
for k = 1:Nit
    
     err_old        = err;
    %%%%%% v-subproblem %%%%%%
    omega = mu./beta;
    
    wgt1 = omega*p./(abs(Dx1)+ 0.00001).^(1-(p)); %IRL1 Weight update;
    wgt2 = omega*p./(abs(Dx2) + 0.00001).^(1-(p));% IRL1 Weight update
    
    v_1 =shrink(Dx1,wgt1); 
    v_2 =shrink(Dx2,wgt2);
    
    
   
   
    %%%%%% x-subproblem %%%%%%
    x_old = x;
    rhs   = beta*Dt(v_1,v_2)+ imfilter(f,K,'circular');%+ rho*v ;
    
   
    
    x     = fft2(rhs)./lhs;
    x     = real(ifft2(x));
    
   
    t      =  (k-1)/(k + 2);
    
    y     = x + t*(x -x_old);
   
    
    [Dx1, Dx2]            = D(y);
    
    
   
    relError_old    = relError(k);
    err            = norm(x - x_old,'fro')/norm(x, 'fro');
    relError(k)    = err;
    

    
    if relError(k) < tol
         break;
    end
    
       
end
tg = toc(tg);

    out.sol                 = x;
    out.relativeError       = relError(1:k);
    out.psnrGain            = psnrGain(1:k);
    out.ssimGain            = ssimGain(1:k);
    out.obj                 = objVal(1:k);
    out.cpuTime             = tg;

end

function [D,Dt] = defDDt()
D  = @(U) ForwardDiff(U);
Dt = @(X,Y) Dive(X,Y);
end

function [Dux,Duy] = ForwardDiff(U)
 Dux = [diff(U,1,2), U(:,1,:) - U(:,end,:)];
 Duy = [diff(U,1,1); U(1,:,:) - U(end,:,:)];
end

function DtXY = Dive(X,Y)
  % Transpose of the forward finite difference operator
  % is the divergence fo the forward finite difference operator
  DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
  DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];   
end


function z = shrink(x,r)
z = sign(x).*max(abs(x)-r,0);
end

