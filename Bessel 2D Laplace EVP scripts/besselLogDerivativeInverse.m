function inverseValue = besselLogDerivativeInverse(nIndex,value,intervalIndex,besselLogDerivativeFncHandle,besselZeroDataTable)

        fAbs = @(y) abs(besselLogDerivativeFncHandle(y)-value);

        if intervalIndex == 1

            inverseValue = fminbnd(fAbs, 0, besselZeroDataTable(nIndex,intervalIndex)-0.000001);

        
        else

            inverseValue = fminbnd(fAbs, besselZeroDataTable(nIndex,intervalIndex-1)+0.000001, besselZeroDataTable(nIndex,intervalIndex)-0.000001);

        end

    end