import tensorflow as tf
import time
from matrix_exp_op import matrix_exp
import sys


def compute_graph_optimize(size, nsteps, n_channels, dt, relax):
    start = time.time()
    # return
    dims = 'xy'
    measure_time = True;
    print("\n\n            Tensorflow is computing the graph:\n\n")

    with tf.Session() as sess:

        zero = tf.zeros([1], dtype=tf.float64)
        zero_vect = tf.zeros([size, 1], dtype=tf.complex128)
        one = tf.ones([1], dtype=tf.float64)
        delta_t = tf.complex(zero, -tf.constant(dt, shape=[1], dtype=tf.float64))
        L0 = tf.placeholder(tf.complex128, shape=[size, size], name='L0')
        if relax:
            R = tf.multiply(tf.complex(zero, one), tf.placeholder(tf.complex128, shape=[size, size], name='R'))  #
            L = L0 + R
        else:
            L = L0
        L_list = [tf.placeholder(tf.complex128, shape=[size, size], name='L' + dim + str(c)) for c in range(n_channels)
                  for dim in dims]
        psi_list = [None] * (nsteps + 1)
        psi_list[0] = tf.placeholder(tf.complex128, shape=[size, 1], name='psi0')
        target_state = tf.placeholder(tf.complex128, shape=[size, 1], name='target')

        U_list = [tf.Variable(tf.placeholder(tf.float64, shape=[1], name="Assign" + dim + str(c) + "_" + str(i)),
                              name="U" + dim + str(c) + "_" + str(i), dtype=tf.float64) for i in range(nsteps) for c in
                  range(n_channels) for dim in dims]
        init = tf.variables_initializer(tf.global_variables(), name='init_all_vars_op')

        if measure_time:
            print("Time needed for the initialization: ", time.time() - start)
            start = time.time()

        def propagate(my_u_list, n):
            hamilt_list = [tf.multiply(
                tf.complex(my_u_list[c * 2 + counter], zero, name="Complex_U" + dim + str(c) + "_" + str(n)),
                L_list[c * 2 + counter], name='Mult' + dim + str(c) + "_" + str(n)) for counter, dim in enumerate(dims)
                for c in range(n_channels)]  # tf.multiply
            hamilt_list.append(L)
            tot_hamil = tf.add_n(hamilt_list, name='Tot_Hamiltonian')
            return matrix_exp(tot_hamil, delta_t=dt, name='Exp' + "_" + str(n))

        for n in range(0, nsteps):
            propagation = propagate(U_list[n * n_channels * 2:((n + 1) * n_channels * 2 + 1)], n)
            psi_list[n + 1] = tf.matmul(propagation, psi_list[n], name='psi' + str(n + 1))

        if measure_time:
            print("Time needed for creating the forward nodes: ", time.time() - start)
            start = time.time()

        final_state = psi_list[nsteps]
        trace = tf.matmul(tf.transpose(final_state, conjugate=True), target_state, name="trace")

        # spin-scenario: write your own fidelity function here.

        #fidelity = tf.real(-trace, name="trace_real")
        fidelity = tf.add(zero, -tf.real(trace), name="trace_real")
        # fidelity = tf.add(zero, -tf.real(trace), name="trace_real")
        #fidelity = tf.add(tf.reduce_sum(tf.nn.l2_normalize(U_list)), -tf.real(trace), name="trace_real")
        deriv = tf.gradients(fidelity, U_list, name="Cecilia_gradient")

        if measure_time:
            print("Time needed for creating the gradient nodes:  ", time.time() - start)

        tf.train.write_graph(sess.graph_def, '../', 'opt_graph.pb', as_text=False)

        # the following lines create the tensorboard log files necessary to view the graph in tensorflow
        merged = tf.summary.merge_all()
        train_writer = tf.summary.FileWriter('tboard/' + '/train', sess.graph)

    return "\n\n            Tensorflow graph succesfully created\n "


def compute_graph_evaluation(size, nsteps, n_channels, dt, relax):  # , relax = False, target
    # return
    dims = 'xy'
    measure_time = True;

    print("\n\n            Tensorflow is computing the graph:\n\n")

    with tf.Session() as sess:

        zero = tf.zeros([1], dtype=tf.float64)
        zero_vect = tf.zeros([size, 1], dtype=tf.complex128)
        one = tf.ones([1], dtype=tf.float64)
        delta_t = tf.complex(zero, -tf.constant(dt, shape=[1], dtype=tf.float64))
        L0 = tf.placeholder(tf.complex128, shape=[size, size], name='L0')
        if relax:
            R = tf.multiply(tf.complex(zero, one), tf.placeholder(tf.complex128, shape=[size, size], name='R'))  #
            L = L0 + R
        else:
            L = L0
        L_list = [tf.placeholder(tf.complex128, shape=[size, size], name='L' + dim + str(c)) for c in range(n_channels)
                  for dim in dims]
        psi_list = [None] * (nsteps + 1)
        psi_list[0] = tf.placeholder(tf.complex128, shape=[size, 1], name='psi0')
        target_state = tf.placeholder(tf.complex128, shape=[size, 1], name='target')
        # here we can make U_list a list of Placeholders! much better!
        U_list = [tf.Variable(tf.placeholder(tf.float64, shape=[1], name="Assign" + dim + str(c) + "_" + str(i)),
                              name="U" + dim + str(c) + "_" + str(i), dtype=tf.float64) for i in range(nsteps) for c in
                  range(n_channels) for dim in dims]
        init = tf.variables_initializer(tf.global_variables(), name='init_all_vars_op')

        if measure_time:
            print("Time needed for the initialization: ", time.time() - start)
            start = time.time()

        def body_inside(rho, next_term, k, iter, my_L):
            help = tf.add(k, one)
            help2 = (delta_t / tf.to_complex128(tf.complex(help * iter, zero)))
            temp = tf.matmul(help2 * my_L, next_term)
            return [rho + temp, temp, help, iter, my_L]

        def condition_inside(rho, next_term, k, iter, my_L):
            term = next_term
            temp = tf.maximum(tf.to_double(sys.float_info.min), tf.reduce_max(tf.sqrt(tf.abs(term))))
            return temp > sys.float_info.epsilon

        def body_outside(rho, next_term, H_tot, iter, ii):
            result = tf.while_loop(condition_inside, body_inside, [zero_vect, next_term, zero, iter, H_tot])
            return [result[0], result[1], H_tot, iter, ii + 1]

        def condition_outside(rho, next_term, H_tot, iter, ii):
            return tf.reshape(tf.greater(iter, ii), [])

        def step(rho0, H_tot, n):
            with tf.name_scope("my_matrix_exp") as scope:
                # rho=expm(-1i*L*dt)*rho;

                # Get the norm of the vector
                norm_rho = tf.norm(rho0)
                # Scale the vector
                rho = rho0 / norm_rho

                # Get the norm of the i*L*dt matrix
                norm_mat = tf.norm(H_tot, ord=1) * tf.abs(delta_t)
                # Determine the number of time steps
                iter = tf.ceil(norm_mat / 5)

                # Run the Krylov procedure
                # result = tf.while_loop(condition_outside, body_outside,[zero_vect, rho, H_tot, iter, zero])
                result = tf.while_loop(condition_outside, body_outside, [zero_vect, rho, H_tot, iter, zero])

            # Scale the vector back
            return tf.multiply(result[0], norm_rho, name='psi' + str(n + 1))

        def get_matvecexp(my_u_list, n):
            # matrix vector exponential
            hamilt_list = [tf.multiply(
                tf.complex(my_u_list[c * 2 + counter], zero, name="Complex_U" + dim + str(c) + "_" + str(n)),
                L_list[c * 2 + counter], name='Mult' + dim + str(c) + "_" + str(n)) for counter, dim in enumerate(dims)
                for c in range(n_channels)]  # tf.multiply
            hamilt_list.append(L)
            tot_hamil = tf.add_n(hamilt_list, name='Tot_Hamiltonian')

            return step(psi_list[n], tot_hamil, n)

        for n in range(0, nsteps):
            psi_list[n + 1] = get_matvecexp(U_list[n * n_channels * 2:((n + 1) * n_channels * 2 + 1)], n)

        if measure_time:
            print("Time needed for creating the forward nodes: ", time.time() - start)
            start = time.time()

        final_state = psi_list[nsteps]
        trace = tf.matmul(tf.transpose(final_state, conjugate=True), target_state, name="trace")
        fidelity = tf.subtract(tf.constant(2.0, dtype=tf.float64), tf.real(trace, name="trace_real"),
                               name="fidelity")  # fidelity for tensorflow

        tf.train.write_graph(sess.graph_def, '../build', 'Mygraph.pb', as_text=False)

        # the following lines create the tensorboard log files necessary to view the graph in tensorflow
        merged = tf.summary.merge_all()
        train_writer = tf.summary.FileWriter('tboard/' + '/train', sess.graph)

    return "\n\n            Tensorflow graph succesfully created\n "
