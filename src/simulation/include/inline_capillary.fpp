#:def compute_capillary_stress_tensor()

    Omega(1, 1) = -sigma_coeff*(w2*w2 + w3*w3)/normW

    Omega(2, 1) = sigma_coeff*w1*w2/normW
    Omega(1, 2) = Omega(2, 1)

    Omega(2, 2) = -sigma_coeff*(w1*w1 + w3*w3)/normW

    if (p > 0) then
        #:if not MFC_CASE_OPTIMIZATION or num_dims > 2
            Omega(3, 1) = sigma_coeff*w1*w3/normW
            Omega(1, 3) = Omega(3, 1)

            Omega(3, 2) = sigma_coeff*w2*w3/normW
            Omega(2, 3) = Omega(3, 2)

            Omega(3, 3) = -sigma_coeff*(w1*w1 + w2*w2)/normW
        #:endif

    end if

#:enddef compute_capillary_stress_tensor
