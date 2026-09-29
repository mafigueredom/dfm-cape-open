using System;
using System.Collections.Generic;
using System.Globalization;
using System.Windows.Forms;

namespace PhD.DfmMethanation
{
    /// <summary>
    /// Modal editor for CAPE-OPEN input parameters. COFE calls ICapeUtilities.Edit
    /// from the Parameters page; E_NOTIMPL surfaces as "Not implemented".
    /// </summary>
    sealed class ParameterEditorForm : Form
    {
        readonly IList<CapeParameter> _inputs;
        readonly DataGridView _grid = new DataGridView();

        public ParameterEditorForm(IList<CapeParameter> inputs)
        {
            _inputs = inputs;
            Text = "DFM Methanation Reactor";
            StartPosition = FormStartPosition.CenterScreen;
            MinimizeBox = false;
            MaximizeBox = true;
            ShowInTaskbar = false;
            Width = 720;
            Height = 520;
            MinimumSize = new System.Drawing.Size(480, 320);

            _grid.Dock = DockStyle.Fill;
            _grid.AllowUserToAddRows = false;
            _grid.AllowUserToDeleteRows = false;
            _grid.AllowUserToResizeRows = false;
            _grid.RowHeadersVisible = false;
            _grid.SelectionMode = DataGridViewSelectionMode.FullRowSelect;
            _grid.MultiSelect = false;
            _grid.AutoSizeColumnsMode = DataGridViewAutoSizeColumnsMode.Fill;
            _grid.Columns.Add(ReadOnlyColumn("Parameter", 28));
            _grid.Columns.Add(new DataGridViewTextBoxColumn
            {
                Name = "Value",
                HeaderText = "Value",
                FillWeight = 24
            });
            _grid.Columns.Add(ReadOnlyColumn("Description", 48));

            foreach (var p in _inputs)
            {
                var note = p.ComponentDescription ?? "";
                if (p is CapeOptionParameter opt && opt.RestrictedToList && opt.Choices.Length > 0)
                    note = note + " [" + string.Join(", ", opt.Choices) + "]";
                int row = _grid.Rows.Add(
                    p.ComponentName,
                    Convert.ToString(p.value, CultureInfo.InvariantCulture),
                    note);
                _grid.Rows[row].Tag = p;
            }

            var buttons = new FlowLayoutPanel
            {
                Dock = DockStyle.Bottom,
                FlowDirection = FlowDirection.RightToLeft,
                Height = 40,
                Padding = new Padding(8)
            };
            var ok = new Button { Text = "OK", DialogResult = DialogResult.None, Width = 88 };
            var cancel = new Button { Text = "Cancel", DialogResult = DialogResult.Cancel, Width = 88 };
            ok.Click += OnOk;
            buttons.Controls.Add(cancel);
            buttons.Controls.Add(ok);
            AcceptButton = ok;
            CancelButton = cancel;

            Controls.Add(_grid);
            Controls.Add(buttons);
        }

        static DataGridViewColumn ReadOnlyColumn(string header, float weight)
        {
            return new DataGridViewTextBoxColumn
            {
                HeaderText = header,
                ReadOnly = true,
                FillWeight = weight
            };
        }

        void OnOk(object sender, EventArgs e)
        {
            _grid.EndEdit();
            var pending = new List<KeyValuePair<CapeParameter, object>>();
            foreach (DataGridViewRow row in _grid.Rows)
            {
                var p = (CapeParameter)row.Tag;
                var text = Convert.ToString(row.Cells[1].Value) ?? "";
                if (!TryParse(p, text.Trim(), out var parsed, out var error))
                {
                    MessageBox.Show(this, error, Text, MessageBoxButtons.OK, MessageBoxIcon.Warning);
                    return;
                }
                pending.Add(new KeyValuePair<CapeParameter, object>(p, parsed));
            }
            foreach (var kv in pending)
                kv.Key.value = kv.Value;
            DialogResult = DialogResult.OK;
            Close();
        }

        static bool TryParse(CapeParameter p, string text, out object parsed, out string error)
        {
            parsed = null;
            error = null;
            if (p is CapeBooleanParameter)
            {
                if (text == "1" || text.Equals("true", StringComparison.OrdinalIgnoreCase))
                {
                    parsed = true;
                    return true;
                }
                if (text == "0" || text.Equals("false", StringComparison.OrdinalIgnoreCase))
                {
                    parsed = false;
                    return true;
                }
                error = p.ComponentName + " must be true or false.";
                return false;
            }
            if (p is CapeIntegerParameter n)
            {
                if (!int.TryParse(text, NumberStyles.Integer, CultureInfo.InvariantCulture, out var i) &&
                    !int.TryParse(text, NumberStyles.Integer, CultureInfo.CurrentCulture, out i))
                {
                    error = p.ComponentName + " must be an integer.";
                    return false;
                }
                string message = "";
                if (!n.Validate(i, ref message))
                {
                    error = message;
                    return false;
                }
                parsed = i;
                return true;
            }
            if (p is CapeOptionParameter opt)
            {
                if (opt.RestrictedToList)
                {
                    string message = "";
                    if (!opt.Validate(text, ref message))
                    {
                        error = message;
                        return false;
                    }
                }
                parsed = text;
                return true;
            }
            if (!double.TryParse(text, NumberStyles.Float, CultureInfo.CurrentCulture, out var d) &&
                !double.TryParse(text, NumberStyles.Float, CultureInfo.InvariantCulture, out d))
            {
                error = p.ComponentName + " must be a number.";
                return false;
            }
            parsed = d;
            return true;
        }
    }
}
