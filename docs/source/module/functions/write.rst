:orphan:

**************
write
**************

In principle, all kinds of data in Spin-Scenario can be save into files by ``write`` function.

.. code-block:: C++

  // This is for sequence blocks.
  void write(string file, sol::variadic_args va, const seq_block & /*sb*/);

  // This is for kinds of matrix or vector in Eigen format.
  void write(string file, sol::variadic_args va, const mat & /*m*/);
  void write(string file, sol::variadic_args va, const cx_mat & /*m*/);
  void write(string file, sol::variadic_args va, const vec & /*v*/);
  void write(string file, sol::variadic_args va, const cx_vec & /*v*/);
  void write(string file, sol::variadic_args va, const sp_mat & /*m*/);
  void write(string file, sol::variadic_args va, const sp_cx_mat & /*m*/);
  void write(string file, sol::variadic_args va, const sp_vec & /*v*/);
  void write(string file, sol::variadic_args va, const sp_cx_vec & /*v*/);

Examples of usage can be found in :doc:`table2vec<table2vec>`, :doc:`table2mat<table2mat>`, :doc:`hardRF<hardRF>`, :doc:`shapedRF<shapedRF>`.