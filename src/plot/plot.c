#include "plot.h"

static int plot_python_initialized = 0;

static int err_py(PyObject *obj)
{
  if (!obj)
  {
    PyErr_Print();
    return -1;
  }

  return 0;
}

static int plot_python_init(void)
{
  if (plot_python_initialized)
    return 0;

  Py_Initialize();

  if (!Py_IsInitialized())
    return -1;

  if (_import_array() < 0)
  {
    PyErr_Print();
    Py_Finalize();
    return -1;
  }

  PyRun_SimpleString(
      "import sys, os; "
      "sys.path.insert(0, os.path.expanduser('./src/plot'))"
  );

  plot_python_initialized = 1;

  return 0;
}

static void plot_python_finalize(void)
{
  if (!plot_python_initialized)
    return;

  Py_Finalize();

  plot_python_initialized = 0;
}

static int plot_python_call(
    const char *func_name,
    PyObject *args
)
{
  int status = -1;

  PyObject *name   = NULL;
  PyObject *module = NULL;
  PyObject *func   = NULL;
  PyObject *result = NULL;

  if (plot_python_init() != 0) goto cleanup;

  name = PyUnicode_FromString(PLOT_MODULE);

  if (err_py(name) != 0)
    goto cleanup;

  module = PyImport_Import(name);

  if (err_py(module) != 0) goto cleanup;

  func = PyObject_GetAttrString(module, func_name);

  if (err_py(func) != 0) goto cleanup;

  if (!PyCallable_Check(func))
  {
    PyErr_Format(
        PyExc_TypeError,
        "%s is not callable",
        func_name
    );

    PyErr_Print();
    goto cleanup;
  }

  result = PyObject_CallObject(func, args);

  if (err_py(result) != 0) goto cleanup;

  status = 0;

cleanup:
  Py_XDECREF(result);
  Py_XDECREF(func);
  Py_XDECREF(module);
  Py_XDECREF(name);

  return status;
}

int plot_model(model_t* model)
{
  int status = -1;

  PyObject *py_model = NULL;
  PyObject *args     = NULL;

  if (plot_python_init() != 0) goto cleanup;

  npy_intp dims[2] = {model->nz, model->nx};

  py_model = PyArray_SimpleNewFromData(
    2,
    dims,
    NPY_FLOAT32,
    (void *)model->vp
  );

  if (err_py(py_model) != 0) goto cleanup;

  args = PyTuple_Pack(1, py_model);

  if (err_py(args) != 0) goto cleanup;

  status = plot_python_call("plot_model", args);

cleanup:
  Py_XDECREF(args);
  Py_XDECREF(py_model);

  return status;
}

int plot_seismogram(seismogram_t* seismogram, int offset)
{
  seismogram_t* s = seismogram;

  int status = -1;

  PyObject *py_seismogram = NULL;
  PyObject *dt_obj        = NULL;
  PyObject *offset_obj    = NULL;
  PyObject *args          = NULL;

  if (plot_python_init() != 0) goto cleanup;

  npy_intp dims[2] = {s->nt, s->nrec};

  py_seismogram = PyArray_SimpleNewFromData(
    2,
    dims,
    NPY_FLOAT32,
    (void *)s->seismogram
  );

  if (err_py(py_seismogram) != 0) goto cleanup;

  dt_obj = PyFloat_FromDouble(s->dt);

  if (err_py(dt_obj) != 0) goto cleanup;

  offset_obj = PyLong_FromLong(offset);

  if (err_py(offset_obj) != 0) goto cleanup;

  args = PyTuple_Pack(
      3,
      py_seismogram,
      dt_obj,
      offset_obj
  );

  if (err_py(args) != 0) goto cleanup;

  status = plot_python_call("plot_seismogram", args);

cleanup:
  Py_XDECREF(args);
  Py_XDECREF(offset_obj);
  Py_XDECREF(dt_obj);
  Py_XDECREF(py_seismogram);

  return status;
}

int plot_model_geometry(model_t* model, int dh, geometry_t* geometry)
{
  model_t* m    = model;
  geometry_t* g = geometry;

  int status = -1;

  PyObject *py_model = NULL;
  PyObject *py_recx  = NULL;
  PyObject *py_recz  = NULL;
  PyObject *py_srcx  = NULL;
  PyObject *py_srcz  = NULL;
  PyObject *py_nb    = NULL;
  PyObject *py_dh    = NULL;
  PyObject *args     = NULL;

  if (plot_python_init() != 0) goto cleanup;

  npy_intp model_dims[2] = {m->nzz, m->nxx};

  py_model = PyArray_SimpleNewFromData(
    2,
    model_dims,
    NPY_FLOAT32,
    (void *)m->vp
  );

  if (err_py(py_model) != 0) goto cleanup;

  npy_intp rec_dims[1] = {g->nrec};

  py_recx = PyArray_SimpleNewFromData(
      1,
      rec_dims,
      NPY_FLOAT32,
      (void *)g->rec.x
  );

  if (err_py(py_recx) != 0) goto cleanup;

  py_recz = PyArray_SimpleNewFromData(
      1,
      rec_dims,
      NPY_FLOAT32,
      (void *)g->rec.z
  );

  if (err_py(py_recz) != 0) goto cleanup;

  npy_intp src_dims[1] = {g->nsrc};

  py_srcx = PyArray_SimpleNewFromData(
      1,
      src_dims,
      NPY_FLOAT32,
      (void *)g->src.x
  );

  if (err_py(py_srcx) != 0) goto cleanup;

  py_srcz = PyArray_SimpleNewFromData(
      1,
      src_dims,
      NPY_FLOAT32,
      (void *)g->src.z
  );

  if (err_py(py_srcz) != 0) goto cleanup;

  py_nb = PyLong_FromLong(m->nb);

  if (err_py(py_nb) != 0) goto cleanup;

  py_dh = PyLong_FromLong(dh);

  if (err_py(py_dh) != 0) goto cleanup;

  args = PyTuple_Pack(
      7,
      py_model,
      py_nb,
      py_dh,
      py_recx,
      py_recz,
      py_srcx,
      py_srcz
  );

  if (err_py(args) != 0) goto cleanup;

  status = plot_python_call("plot_model_geometry", args);

cleanup:
  Py_XDECREF(args);

  Py_XDECREF(py_model);
  Py_XDECREF(py_recx);
  Py_XDECREF(py_recz);
  Py_XDECREF(py_srcx);
  Py_XDECREF(py_srcz);

  Py_XDECREF(py_nb);
  Py_XDECREF(py_dh);

  return status;
}

int plot1d(float* arr, int size)
{
  int status = -1;

  PyObject *py_arr = NULL;
  PyObject *args     = NULL;

  if (plot_python_init() != 0) goto cleanup;

  npy_intp dims[1] = {size};

  py_arr = PyArray_SimpleNewFromData(
    2,
    dims,
    NPY_FLOAT32,
    arr
  );

  if (err_py(py_arr) != 0) goto cleanup;

  args = PyTuple_Pack(1, py_arr);

  if (err_py(args) != 0) goto cleanup;

  status = plot_python_call("plot1d", args);

cleanup:
  Py_XDECREF(args);
  Py_XDECREF(py_arr);

  return status;
}

int plot2d(float* arr, int row, int col)
{
  int status = -1;

  PyObject *py_arr = NULL;
  PyObject *args     = NULL;

  if (plot_python_init() != 0) goto cleanup;

  npy_intp dims[2] = {row, col};

  py_arr = PyArray_SimpleNewFromData(
    2,
    dims,
    NPY_FLOAT32,
    arr
  );

  if (err_py(py_arr) != 0) goto cleanup;

  args = PyTuple_Pack(1, py_arr);

  if (err_py(args) != 0) goto cleanup;

  status = plot_python_call("plot2d", args);

cleanup:
  Py_XDECREF(args);
  Py_XDECREF(py_arr);

  return status;
}
