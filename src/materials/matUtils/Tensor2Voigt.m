function eVoigt = Tensor2Voigt(eTensor,dim)

switch dim

    case 2

        eVoigt         = zeros(dim,1);
        eVoigt(1,1)    = eTensor(1,1);
        eVoigt(2,1)    = eTensor(2,2);
        eVoigt(3,1)    = 0.5 * ( eTensor(1,2) + eTensor(2,1) );

    case 3

        error('Not yet implemented!')

    otherwise
        error('Dimension and array mismatch!')
end

end