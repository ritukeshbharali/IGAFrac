function eTensor = Voigt2Tensor(eVoigt,dim,trueStrain)

switch dim

    case 2

        eTensor       = zeros(dim);
        eTensor(1,1)  = eVoigt(1);
        eTensor(2,2)  = eVoigt(2);

        if trueStrain == true
            eVoigt(3) = eVoigt(3)/2;
        end

        eTensor(1,2)  = eVoigt(3);
        eTensor(2,1)  = eVoigt(3);

    case 3

        error('Not yet implemented!')

    otherwise
        error('Dimension and array mismatch!')
end

end