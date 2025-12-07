function [JmZeroK,JofJ]=BessDerivZerosBisect2(nVector,mVector,tol)
%BessDerivZerosBisect2 Zeros of the first derivative of Bessel function of the first kind.
% function [jprimemk,JofJ]=BessDerivZerosBisect2(MM,KK,tol)
% Calculates the zeros of the first derivatives of Bessel function J'm.
%  MM = a vector of the Bessel function orders (>= 0)
%  KK = a vector of the zero indices (>= 1)
%  tol = tolerance on the value of the derivative, not on the value of the zero
%        (May hang for tol < 1e-12)
%
%Use the bisection algorithm, because it gives the desired roots.
%Estimates for bracketing intervals are taken from the
%asymptotic forms for the roots, in Abramowitz and Stegun.
% Following tradition, only positive zeros are given.  
% So the first zero for J0 is 3.83...
%
% jprimemk = the roots
% JofJ     = Bessel function evaluated there.
%
% Routine written by Larry Forbes, University of Tasmania.
% Modified 2010/06/23  Carey Smith
% Modified 2010/07/09  Carey Smith--Changed the initial guess for large values of m
%                      Carey Smith--Corrected the initial guess for n=17, m=5
% Modified 2011/4/27   Vincent: "It tends to crash for high values of m and k > 1.
%                           (cf. for instance m=44 and k=4)
%                            I fixed it, though, by replacing lines 103 to 110."

if nargin == 2
    tol = 1.e-12;
end
thetaIndicies = length(nVector);
rIndicies  = length(mVector);
JmZeroK = zeros(thetaIndicies,rIndicies);
JofJ     = JmZeroK;

% Use a table for accuracy and speed. Is is a combination of data from:
%   http://wwwal.kuicr.kyoto-u.ac.jp/www/accelerator/a4/besselroot.htmlx
%   Abramowitz and Stegun (p.411)
%   and this routine
% It is transposed from A&S:
% Rows are values of m (starting at 0); columns are values of k (starting at 1)
%    k=1              k=2              k=3              k=4              k=5              k=6      k=7      k=8      k=9      k=10     k=11     k=12     k=13     k=14     k=15     k=16     k=17     k=18     k=19     k=20
BesselDerivativeZerosT = [...
    3.83170597020751 7.01558666981561 10.1734681350627 13.3236919363142 16.4706300508776 19.6158585105 22.7600843806 25.9036720876 29.0468285349 32.1896799110 35.3323075501 38.4747662348 41.6170942128 44.7593189977 47.9014608872 51.0435351836 54.1855536411 57.3275254379 60.4694578453 63.61136 % m=0
    1.84118378134065 5.33144277352503 8.53631636634628 11.7060049025920 14.8635886339090 18.01553 21.16437 24.31133 27.45705 30.60192 33.74618 36.88999 40.03344 43.17663 46.31960 49.46239 52.60504 55.74757 58.89000 62.03235 % m=1
    3.05423692822714 6.70613319415845 9.96946782308759 13.1703708560161 16.3475223183217 19.51291 22.67158 25.82604 28.97767 32.12733 35.27554 38.42265 41.56893 44.71455 47.85964 51.00430 54.14860 57.29260 60.43635 63.57989 % m=2
    4.20118894121052 8.01523659837595 11.3459243107430 14.5858482861670 17.7887478660664 20.97248 24.14490 27.31006 30.47027 33.62695 36.78102 39.93311 43.08365 46.23297 49.38130 52.52882 55.67567 58.82195 61.96775 65.11315 % m=3
    5.31755312608399 9.28239628524161 12.6819084426388 15.9641070377315 19.1960288000489 22.40103 25.58976 28.76784 31.93854 35.10392 38.26532 41.42367 44.57962 47.73367 50.88616 54.03737 57.18752 60.33677 63.48526 66.63309 % m=4
    6.41561637570024 10.5198608737723 13.9871886301403 17.3128424878846 20.5755145213868 23.80358 27.01031 30.20285 33.38544 36.56078 39.73064 42.89627 46.05857 49.21817 52.37559 55.53120 58.68528 61.83809 64.98980 68.14057 % m=5
    7.50126614468414 11.7349359530427 15.2681814610978 18.6374430096662 21.9317150178022 25.18393 28.40978 31.61788 34.81339 37.99964 41.17885 44.35258 47.52196 50.68782 53.85079 57.01138 60.16995 63.32681 66.48221 69.63635 % m=6
    8.57783648971407 12.9323862370895 16.5293658843669 19.9418533665273 23.2680529264575 26.54503 29.79075 33.01518 36.22438 39.42227 42.61152 45.79400 48.97107 52.14375 55.31282 58.47887 61.64239 64.80374 67.96324 71.12113 % m=7
    9.64742165199721 14.1155189078946 17.7740123669152 21.2290626228531 24.5871974863176 27.88927 31.15533 34.39663 37.62008 40.83018 44.03001 47.22176 50.40702 53.58700 56.76260 59.93454 63.10340 66.26961 69.43356 72.59554 % m=8
    10.7114339706999 15.2867376673329 19.0045935379460 22.5013987267772 25.8912772768391 29.21856 32.50525 35.76379 39.00190 42.22464 45.43548 48.63692 51.83078 55.01844 58.20095 61.37915 64.55368 67.72509 70.89378 74.06014 % m=9
    11.7708766749555 16.4478527484865 20.2230314126817 23.7607158603274 27.1820215271905 30.53451 33.84197 37.11800 40.37107 43.60677 46.82896 50.04043 53.24322 56.43889 59.62863 62.81338 65.99389 69.17075 72.34447 75.51545 % m=10
    ];

%besselFncZeroStartingTable = [...
%    2.4048255577 5.5200781103 8.6537279129 11.7915344391 14.9309177086 18.0710639679 21.2116366299 24.3524715308 27.4934791320 30.634606484 33.7758202136 36.9170983537 40.0584257646 43.1997917132 46.3411883717 49.4826098974 52.6240518411 55.7655107550 58.9069839261 62.0484691902 % m=0
%    3.83170597020751 7.01558666981561 10.1734681350627 13.3236919363142 16.4706300508776 19.6158585105 22.7600843806 25.9036720876 29.0468285340 32.1896799110 35.3323075501 38.4747662348 41.6170942128 44.7593189977 47.9014608872 51.0435351836 54.1855536411 57.3275254379 60.4694578453 63.6113566985 % m=1
%    5.135622 8.417244 11.619841 14.795952 17.959819 21.116997 24.270112 27.420574 30.569205 33.716520 36.86286 40.00845 43.15345 46.29800 49.44216 51.00430 54.14860 57.29260 60.43635 63.57989 % m=2
%    6.38016 9.76102 13.01520 16.22347 19.40942 22.58273 25.74817 28.90835 32.06485 35.21867 38.37047 41.52072 44.66974 47.81779 50.96503 54.11162 57.25765 60.40322 63.54840 66.69324 %m=3
%    ];

besselFncZeroStartingTable = [...
    2.4048255577 5.5200781103 8.6537279129 11.7915344391 14.9309177086 18.0710639679 21.2116366299 24.3524715308 27.4934791320 30.634606484 33.7758202136 36.9170983537 40.0584257646 43.1997917132 46.3411883717 49.4826098974 52.6240518411 55.7655107550 58.9069839261 62.0484691902 % m=0
    ];
  
% n = 17, abs(m)=5:  56.68528
  
[dataTableMaxRow, dataTableMaxCol] = size(besselFncZeroStartingTable);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nIndex=1:thetaIndicies

    nValue = nVector(nIndex);

    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    for mIndex=1:rIndicies

        mValue = mVector(mIndex);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if((nIndex <= dataTableMaxRow) && (mIndex <= dataTableMaxCol))
            asymptroot = besselFncZeroStartingTable(nIndex,mIndex);
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            if(mIndex <= 5)
                dx1 = 1e-8;  % small, because these table values should be very close
            else
                dx1 = 3e-5;  % small, because these table values should be very close
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            ALeft = asymptroot-dx1;
            ARight= asymptroot+dx1;
        else
            %if(nIndex == 0)
            %    mIndex = mIndex+1;
            %end

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            if mIndex == 1
                oneterm = nValue + (1.8557571*(nValue^(1/3))) + (1.033150 * (nValue^((-1)*(1/3))));
                oneterm = oneterm - (0.00397*(nValue^(-1))) - 0.0908 * (nValue^((-1)*(5/3)));
                oneterm = oneterm + 0.043 * (nValue^((-1)*(7/3)));
                asymptroot = oneterm;
            else
                %k = k+1;
                alphaApr = (mValue + ((1/2)*nValue) - (1/4))*pi;
                mu=4*(nValue^2);
                asymptroot=alphaApr-((mu-1)/(8*alphaApr));
                asymptroot=asymptroot-((4*((mu-1)*(7*mu-31)))/(3*((8*alphaApr)^3)));
                asymptroot = asymptroot - ( 32 * ( ((mu-1) * (83*(mu^2)-982*mu+3779)) / (15*(8*alphaApr)^5) ) );
                %asymptroot = asymptroot - 64*((mu-1)*(6949*(mu^3)-153855*(mu^2)+1585743*mu-6277237))/(105*((8*mu)^7));
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            ALeft = asymptroot - 1;
            ARight = asymptroot + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        JL=besselj(nValue,ALeft);
        JR=besselj(nValue,ARight);
        
        %if JprimeL*JprimeR > 0
        %disp(['BessDerivZerosBisect1 error for m=',int2str(m),', k=',int2str(k)])
        %disp(['   ARight=',num2str(ARight),', asymptroot=',num2str(asymptroot),', ALeft=',num2str(ALeft)])
        %warning('BessDerivZerosBisect1:Interval_Error',['Original interval does not bracket root, JprimeL=',num2str(JprimeL),', JprimeR=',num2str(JprimeR)]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while JL*JR > 0

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            if(JL > 0)
                ALeft = ALeft - (ARight - ALeft)/2;
                JL=besselj(nValue,ALeft);
            else
                ARight = ARight + (ARight - ALeft)/2;
                JR=besselj(nValue,ARight);
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Amiddle=(ALeft+ARight)/2;
        JM=besselj(nValue,Amiddle);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        while abs(JM) > tol

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            if JL*JM < 0
                ARight=Amiddle;
            else
                ALeft=Amiddle;
                JL=besselj(nValue,ALeft);
            end
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

            Amiddle=(ALeft+ARight)/2;
            JM=besselj(nValue,Amiddle);

        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        JmZeroK(nIndex,mIndex) = Amiddle;
        JofJ(nIndex,mIndex) = JM;

    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end



%tic
%(besselj(2,8)-besselj(4,8))/2;
%toc

%tic
%dx = 0.0001;
%(besselj(3,8+dx) - besselj(3,max(8-dx,0)))/(2*dx)
%toc

