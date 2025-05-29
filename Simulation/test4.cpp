void propagator(double t_max = 0, double dt = 0.001)
    {
        double t = 0.0;
        int steps = static_cast<int>(t_max / dt);

        // Constants for field
        double qmdt2 = (q / m) * (dt / 2.0);

        for (int i = 0; i < steps; i++)
        {
            // Half acceleration from E
            double vx_minus = vx + qmdt2 * Ex;
            double vy_minus = vy + qmdt2 * Ey;
            double vz_minus = vz + qmdt2 * Ez;

            // t vector
            double tx = qmdt2 * Bx;
            double ty = qmdt2 * By;
            double tz = qmdt2 * Bz;

            // t^2
            double t_mag2 = tx * tx + ty * ty + tz * tz;

            // s vector
            double sx = 2 * tx / (1 + t_mag2);
            double sy = 2 * ty / (1 + t_mag2);
            double sz = 2 * tz / (1 + t_mag2);

            // v' = v_minus + v_minus x t
            double vpx = vx_minus + (vy_minus * tz - vz_minus * ty);
            double vpy = vy_minus + (vz_minus * tx - vx_minus * tz);
            double vpz = vz_minus + (vx_minus * ty - vy_minus * tx);

            // v_plus = v_minus + v' x s
            double vx_plus = vx_minus + (vpy * sz - vpz * sy);
            double vy_plus = vy_minus + (vpz * sx - vpx * sz);
            double vz_plus = vz_minus + (vpx * sy - vpy * sx);

            // Final velocity (after second half E kick)
            vx = vx_plus + qmdt2 * Ex;
            vy = vy_plus + qmdt2 * Ey;
            vz = vz_plus + qmdt2 * Ez;

            // Position update using full-step velocity
            x += vx * dt;
            y += vy * dt;
            z += vz * dt;

            // Save position
            posx.push_back(x);
            posy.push_back(y);
            posz.push_back(z);
        }

        // Final velocity and energy
        v = sqrt(vx * vx + vy * vy + vz * vz);
        energy = 0.5 * m * v * v;
    }