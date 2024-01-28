function dist = eucl(x, p)

    diff = x-p;

    diff = diff.^2;
    ss = sum(diff, 1);

    dist= sqrt(ss);

end