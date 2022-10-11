{{ name }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :no-members:
   :no-inherited-members:
   :no-special-members:

   {% block methods %}
   {% if methods %}
   .. HACK -- we don't want this to appear in the output, but autosummary should still generate the pages.
      .. autosummary::
         :toctree:
         {% for item in all_methods %}
            {%- if not item.startswith('_') or item in ['__call__'] %}
            {{ name }}.{{ item }}
            {%- endif -%}
         {%- endfor %}
   {% endif %}
   {% endblock %}
