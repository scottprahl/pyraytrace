.. |pypi| image:: https://img.shields.io/pypi/v/pyraytrace?color=68CA66
   :target: https://pypi.org/project/pyraytrace/
   :alt: PyPI

.. |github| image:: https://img.shields.io/github/v/tag/scottprahl/pyraytrace?label=github&color=68CA66
   :target: https://github.com/scottprahl/pyraytrace
   :alt: GitHub

.. |license| image:: https://img.shields.io/github/license/scottprahl/pyraytrace?color=68CA66
   :target: https://github.com/scottprahl/pyraytrace/blob/main/LICENSE.txt
   :alt: License

.. |test| image:: https://github.com/scottprahl/pyraytrace/actions/workflows/test.yaml/badge.svg
   :target: https://github.com/scottprahl/pyraytrace/actions/workflows/test.yaml
   :alt: Tests

.. |docs| image:: https://readthedocs.org/projects/pyraytrace/badge/?version=latest&color=68CA66
   :target: https://pyraytrace.readthedocs.io
   :alt: Docs

.. |downloads| image:: https://img.shields.io/pypi/dm/pyraytrace?color=68CA66
   :target: https://pypi.org/project/pyraytrace/
   :alt: Downloads

.. |lite| image:: https://img.shields.io/badge/try-JupyterLite-68CA66.svg
   :target: https://scottprahl.github.io/pyraytrace/
   :alt: Try JupyterLite

pyraytrace
==========

|pypi| |github| |license| |test| |docs| |downloads|

|lite|

``pyraytrace`` is a very simple start at Python ray tracing for thin lenses and
basic optical elements. It is intentionally lightweight and educational, with a
focus on clear code and notebook-friendly examples instead of exhaustive optical
modeling.

Status
------

* Early-stage project.
* API is still evolving.
* Best suited for learning and simple experiments.

Features
--------

* 2D paraxial ray tracing tools based on ABCD matrices.
* Simple 3D geometric ray/object primitives (planes, spheres, prisms, lenses).
* Example notebook in ``docs/``.
* Sphinx docs and JupyterLite support.

Installation
------------

Install from PyPI::

   pip install pyraytrace

Quick Example
-------------

.. code-block:: python

   import matplotlib.pyplot as plt
   from pyraytrace.pyparaxial import Object, ParaxialSystem, Plane, Ray, ThinLens

   system = ParaxialSystem()
   system.append(Object(z=0, h=4))
   system.append(ThinLens(z=25, f=35, h=6))
   system.append(Plane(z=65, h=8))

   ray = Ray(y=2, nu=0.0)
   ray.aim(z0=0, y0=2, z1=25, y1=0)

   system.draw()
   system.draw_ray_through_system(ray, color="tab:blue")
   plt.xlabel("z")
   plt.ylabel("y")
   plt.show()

Documentation
-------------

* Read the docs: https://pyraytrace.readthedocs.io


Contributing
------------

Issues and pull requests are welcome at:
https://github.com/scottprahl/pyraytrace/issues

License
-------

``pyraytrace`` is licensed under the terms of the MIT license.
