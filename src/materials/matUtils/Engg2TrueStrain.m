function e = Engg2TrueStrain(e,dim)

switch dim

    case 2

        e(3) = e(3)/2;

    case 3

        e(4) = e(4)/2;
        e(5) = e(5)/2;
        e(6) = e(6)/2;

    otherwise
        error('Wrong dimensions')
end

end