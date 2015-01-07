package com.Sts.UI;

import com.Sts.DBTypes.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsTablePanelNew2 extends StsGroupBox implements ListSelectionListener
{
    private JLabel titleLabel = new JLabel();
    private StsObjectTableModel2 tableModel;
    private JTable table;
    private transient Vector listSelectionListeners = null;
    private StsTableCellRenderer cellRenderer = new StsTableCellRenderer(this);
    private DeleteButtonRendererEditor buttonRendererEditor;
    private int nDeleteButtonColumn = -1;
    private byte[] rowHighlights = null;
    private StsButton deleteButton;
    
    private java.util.List selectRowNotifyListeners = new LinkedList();
    private java.util.List removeRowNotifyListeners = new LinkedList();

    public static final byte NOT_HIGHLIGHTED = 0; // White
    public static final byte HIGHLIGHTED = 1; // Light Blue
    public static final byte NOT_EDITABLE = 2; // Grey
    public static final byte SELECTED = 3; // Blue

    public StsTablePanelNew2()
    {
        try
        {
            buildPanel();
        }
        catch(Exception e)
        {}
    }

    public StsTablePanelNew2(boolean showTitle)
    {
        titleLabel.setVisible(showTitle);
        try
        {
            buildPanel();
        }
        catch(Exception e)
        {}
    }

    public StsTablePanelNew2(Object[] rowObjects, String[] columnNames)
    {
        this(new StsObjectTableModel2(rowObjects, columnNames));
    }

    public StsTablePanelNew2(StsObjectTableModel2 tableModel_)
    {
        final StsObjectTableModel2 tableModel = tableModel_;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                initializer(tableModel);
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    private void initializer(StsObjectTableModel2 tableModel)
    {
        try
        {
            setTableModel(tableModel);
            buildPanel();
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void buildPanel()
    {
        table.setShowGrid(true);
        table.setRequestFocusEnabled(false);
        table.getSelectionModel().addListSelectionListener(this);

        JScrollPane scrollPane = new JScrollPane(table);
        gbc.fill = GridBagConstraints.BOTH;
        scrollPane.setPreferredSize(new Dimension(500, 40));
        add(scrollPane);
    }


    public void deleteRow()
    {
        int selectedRow = table.getSelectedRow();
        System.out.println("Deleting row " + selectedRow);
    }

    public void addMouseListener(MouseListener listener)
    {
        table.addMouseListener(listener);
    }

    public JTable getTable()
    {
        return table;
    }

    public void setTitle(String title)
    {
        titleLabel.setText(title);
    }

    public String getTitle()
    {
        return titleLabel.getText();
    }

    public void setTableModel(StsObjectTableModel2 model)
    {
        tableModel = model;
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                constructTable();
            }
        };
        StsToolkit.runWaitOnEventThread(runnable);
    }

    private void constructTable()
    {
        initializeRowHighlights();
        table = new JTable(tableModel);
        buttonRendererEditor = new DeleteButtonRendererEditor(this);
        setTableRenderer();
        if(tableModel.addDeleteButtons)
            addDeleteButtons();
    }

    private void addDeleteButtons()
    {
        this.nDeleteButtonColumn = tableModel.getColumnCount()-1;
        setTableEditor();
    }

    private void initializeRowHighlights()
    {
        rowHighlights = new byte[tableModel.getRowCount()];
    }

    private void setTableRenderer()
    {
        TableColumn tableColumn;
        int nCols = tableModel.getColumnCount();
        for(int i = 0; i < nCols-1; i++)
        {
            tableColumn = table.getColumnModel().getColumn(i);
            tableColumn.setCellRenderer(cellRenderer);
        }
        if(tableModel.addDeleteButtons)
        {
            tableColumn = table.getColumnModel().getColumn(nCols-1);
            tableColumn.setCellRenderer(buttonRendererEditor);
        }
    }

    private void setTableEditor()
    {
        if(tableModel.addDeleteButtons)
        {
            int nCols = tableModel.getColumnCount();
            TableColumn tableColumn = table.getColumnModel().getColumn(nCols-1);
            tableColumn.setCellEditor(buttonRendererEditor);
        }
    }

    public void addColumns(String[] headers)
    {
        tableModel.initializeColumns(headers);
    }

    public void setAutoResizeMode(int val)
    {
        table.setAutoResizeMode(val);
    }

    public void setColumnWidth(int colIdx, int width)
    {
        TableColumn col = table.getColumnModel().getColumn(colIdx);
        col.setPreferredWidth(width);
    }

    public int getNumberOfRows()
    {
        return table.getRowCount();
    }

    public int[] getSelectedColumns()
    {
        return table.getSelectedColumns();
    }

    public int getNumberOfColumns()
    {
        return table.getColumnCount();
    }
    
    public void setSelectAll(boolean selectAll)
    {
        if( !selectAll) return;

        int[] all = new int[tableModel.getRowCount()];
        for (int i = 0; i < all.length; i++)
           all[i] = i;

        setSelectedIndices(all);
    }

    public void setSelectionIndex(int index)
    {
        setSelectedIndices(new int[] {index});
    }

    public void setColumnSelectable(boolean val)
    {
        table.setColumnSelectionAllowed(val);
    }

    public void setRowsSelectable(boolean val)
    {
        table.setRowSelectionAllowed(val);
    }

    public void setRowType(int row, byte value)
    {
        rowHighlights[row] = value;
//        System.out.println("Setting selectedRows[" + row + "]");
    }

    public void setSelectionMode(int selectionMode)
    {
        table.setSelectionMode(selectionMode);
    }

    public byte getRowHighlight(int row)
    {
        try
        {
            return rowHighlights[row];
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
            return (byte)0;
        }
    }

    private byte[] addElementToArray(byte[] arrayIn, byte value)
    {
        byte[] arrayOut = null;
        if(arrayIn == null)
            arrayOut = new byte[1];
        else
        {
            arrayOut = new byte[arrayIn.length + 1];
            System.arraycopy(arrayIn, 0, arrayOut, 0, arrayIn.length);
        }
        arrayOut[arrayOut.length - 1] = value;
        return arrayOut;
    }

    public int[] getSelectedIndices()
    {
        return table.getSelectedRows();
    }

    public void setColumnSelectable()
    {
        table.setColumnSelectionAllowed(false);
    }

    public void setSelectedIndices(int[] indices)
    {
        DefaultListSelectionModel model = (DefaultListSelectionModel)table.getSelectionModel();
        if(indices == null)
            model.clearSelection();
        else
        {
            for(int i = 0; i < indices.length; i++)
                model.addSelectionInterval(indices[i], indices[i]);
        }
    }

    public void addRow(Object row)
    {
        rowHighlights = addElementToArray(rowHighlights, this.NOT_HIGHLIGHTED);
        tableModel.addRow(row);
        repaint();
    }

    public void addRows(Object[] rows)
    {
        tableModel.addRows(rows);
        initializeRowHighlights();
        table.revalidate();
        repaint();
    }

    public void replaceRows(Object[] rows)
    {
        tableModel.replaceRows(rows);
        initializeRowHighlights();
        repaint();
    }

    public void removeAllRows()
    {
        tableModel.removeAllRows();
        rowHighlights = null;
        repaint();
    }

    public void removeRow(Object row)
    {
        int index = tableModel.removeRowObject(row);
        if(index == StsParameters.NO_MATCH) return;
        rowHighlights = removeElementFromArray(rowHighlights, index);
        repaint();
    }

    private byte[] removeElementFromArray(byte[] arrayIn, int index)
    {
        int newLength = arrayIn.length - 1;
        byte[] arrayOut = new byte[arrayIn.length - 1];
        System.arraycopy(arrayIn, 0, arrayOut, 0, index);
        System.arraycopy(arrayIn, index + 1, arrayOut, index, newLength - index);
        return arrayOut;
    }
    public void rowSelectNotify( StsSelectRowNotifyListener listener, int[] indices)
    {
        listener.rowsSelected( indices);
    }

    public synchronized void addSelectRowNotifyListener(StsSelectRowNotifyListener listener)
    {
        if (! selectRowNotifyListeners.contains(listener))
        {
            selectRowNotifyListeners.add(listener);
        }
        setRowsSelectable( selectRowNotifyListeners.size() > 0);
    }

    public synchronized void removeSelectRowNotifyListener(StsSelectRowNotifyListener listener)
    {
        selectRowNotifyListeners.remove(listener);
        setRowsSelectable( selectRowNotifyListeners.size() > 0);
    }

    public synchronized void addListSelectionListener(ListSelectionListener l)
    {
        Vector v = listSelectionListeners == null ? new Vector(2) : (Vector)listSelectionListeners.clone();
        if(!v.contains(l))
        {
            v.addElement(l);
            listSelectionListeners = v;
        }
    }

    public synchronized void removeListSelectionListener(ListSelectionListener l)
    {
        if(listSelectionListeners != null && listSelectionListeners.contains(l))
        {
            Vector v = (Vector)listSelectionListeners.clone();
            v.removeElement(l);
            listSelectionListeners = v;
        }
    }

    protected void fireValueChanged(ListSelectionEvent e)
    {
        if(listSelectionListeners != null)
        {
            Vector listeners = listSelectionListeners;
            int count = listeners.size();
            for(int i = 0; i < count; i++)
                ((ListSelectionListener)listeners.elementAt(i)).valueChanged(e);
        }
    }

/*
    public void valueChanged(ListSelectionEvent e)
    {
        fireValueChanged(e);
    }
*/
    
    public void valueChanged(ListSelectionEvent e)
    {
        table.requestFocus();
        if(!e.getValueIsAdjusting())
        {
            ListSelectionModel lsm = (ListSelectionModel)e.getSource();
            if (lsm.isSelectionEmpty())
            {
                return;
            }
            int[] indexes = table.getSelectedRows();
            Iterator itr = selectRowNotifyListeners.iterator();
            while (itr.hasNext())
            {
                StsSelectRowNotifyListener listener = (StsSelectRowNotifyListener)itr.next();
                rowSelectNotify(listener, indexes);
            }
        }
    }

    public boolean isSelectedRow(int row)
    {
        for(int i = 0; i < table.getSelectedRowCount(); i++)
        {
            if(row == table.getSelectedRows()[i])
                return true;
        }
        return false;
    }

    public void setValueAt(Object value, Object rowObject, String columnName)
    {
        tableModel.setValueAt(value, rowObject, columnName);
        repaint();
    }

    static public void main(String[] args)
    {

        int nObjects = 3;
        StsSeismicVolume[] rowObjects = new StsSeismicVolume[nObjects];
        String[] columnNames = new String[]{"xMin", "xMax", "yMin", "yMax", "statusString"};
        String[] columnTitles = new String[]{"minX", "maxX", "minY", "maxY", "Status"};
        for(int n = 0; n < nObjects; n++)
        {
            rowObjects[n] = new StsSeismicVolume();
            rowObjects[n].initialize(n + 0.1f, n + 0.2f, n + 0.3f, n + 0.4f, n + 0.5f, n + 0.6f);
            rowObjects[n].statusString = StsSeismicBoundingBox.STATUS_OK_STR;
        }
        StsObjectTableModel2 tableModel = new StsObjectTableModel2(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        StsTablePanelNew2 tablePanel = new StsTablePanelNew2(tableModel);
        StsToolkit.createDialog(tablePanel, false);

        StsToolkit.sleep(2000);
        tablePanel.addRows(rowObjects);

        for(int n = 0; n < nObjects; n++)
        {
            StsToolkit.sleep(500);
            tablePanel.setValueAt(StsSeismicBoundingBox.STATUS_GRID_BAD, rowObjects[n], "statusString");
        }
        for(int n = nObjects - 1; n >= 0; n--)
        {
            StsToolkit.sleep(500);
            tablePanel.removeRow(rowObjects[n]);
        }
        for(int n = 0; n < nObjects; n++)
        {
            StsToolkit.sleep(500);
            tablePanel.addRow(rowObjects[n]);
        }
        StsToolkit.sleep(500);
        tablePanel.addRows(rowObjects);
    }

   class DeleteButtonRendererEditor extends StsButton implements TableCellRenderer, TableCellEditor
    {
        public DeleteButtonRendererEditor(StsTablePanelNew2 panel)
        {
            super("X", panel, "deleteRow");
        }

        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
        {
            return this;
        }

        public boolean isCellEditable(EventObject anEvent)
        {
            return true;
        }

        public Object getCellEditorValue()
        {
            return null;
        }

       public boolean shouldSelectCell(EventObject anEvent)
       {
           return true;
       }

        /**
         * Tells the editor to stop editing and accept any partially edited
         * value as the value of the editor.  The editor returns false if
         * editing was not stopped; this is useful for editors that validate
         * and can not accept invalid entries.
         *
         * @return	true if editing was stopped; false otherwise
         */
        public boolean stopCellEditing()
        {
            return true;
        }

        /**
         * Tells the editor to cancel editing and not accept any partially
         * edited value.
         */
        public void cancelCellEditing()
        {
        }

        /**
         * Adds a listener to the list that's notified when the editor
         * stops, or cancels editing.
         *
         * @param	l		the CellEditorListener
         */
        public void addCellEditorListener(CellEditorListener l)
        {

        }

        /**
         * Removes a listener from the list that's notified
         *
         * @param	l		the CellEditorListener
         */
        public void removeCellEditorListener(CellEditorListener l)
        {

        }

        public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column)
        {
            return this;
        }
    }

    class StsTableCellRenderer extends DefaultTableCellRenderer
    {
        StsTablePanelNew2 panel = null;

        Color bg_highlighted = new Color(255, 255, 165);
        Color fg_highlighted = Color.BLACK;
        Color bg_not_highlighted = Color.WHITE;
        Color fg_not_highlighted = Color.BLACK;
        Color bg_not_editable = new Color(220, 220, 220);
        Color fg_not_editable = Color.BLACK;
        Color bg_selected = Color.BLUE;
        Color fg_selected = Color.WHITE;

        public StsTableCellRenderer(StsTablePanelNew2 panel)
        {
            this.panel = panel;
        }

        // method to override - returns cell renderer component
        public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column)
        {
            if(column == nDeleteButtonColumn)
                return deleteButton;

            // let the default renderer prepare the component for us
            Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);

            // now get the current font used for this cell
            Font font = comp.getFont();
            if(panel.isSelectedRow(row))
            {
                comp.setBackground(bg_selected);
                comp.setForeground(fg_selected);
            }
            else if(panel.getRowHighlight(row) == HIGHLIGHTED)
            {
                comp.setBackground(bg_highlighted);
                comp.setForeground(fg_highlighted);
            }
            else if(panel.getRowHighlight(row) == NOT_EDITABLE)
            {
                comp.setBackground(bg_not_editable);
                comp.setForeground(fg_not_editable);
            }
            else if(panel.getRowHighlight(row) == NOT_HIGHLIGHTED)
            {
                comp.setBackground(bg_not_highlighted);
                comp.setForeground(fg_not_highlighted);
            }
            return comp;
        }
    }
}
