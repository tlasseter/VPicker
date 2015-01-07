//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class StsTablePanel extends StsJPanel implements ListSelectionListener
{
   JLabel titleLabel = new JLabel();
   DefaultTableModel tableModel = new DefaultTableModel();
   public JTable table = new JTable(tableModel);
   JScrollPane scrollPane = new JScrollPane(table);
   private transient Vector listSelectionListeners = null;
   DefaultTableCellRenderer cellRender = new DefaultTableCellRenderer();
   byte[] selectedRows = null;

   public static final byte NOT_HIGHLIGHTED = 0; // White
   public static final byte HIGHLIGHTED = 1; // Light Blue
   public static final byte NOT_EDITABLE = 2; // Grey
   public static final byte SELECTED = 3; // Blue

   public StsTablePanel()
   {
      try
      {
         buildPanel();
      }
      catch (Exception e)
      {}
   }

   public StsTablePanel(boolean showTitle)
   {
      titleLabel.setVisible(showTitle);
      try
      {
         buildPanel();
      }
      catch (Exception e)
      {}
   }

   public StsTablePanel(int nrows, Object[] colNames)
   {
        this(new DefaultTableModel(colNames, nrows));
   }

    public StsTablePanel(TableModel tableModel)
    {
       try
       {
          setTableModel(tableModel);
          table.setModel(tableModel);
          buildPanel();
       }
       catch (Exception e)
       {
          e.printStackTrace();
       }
    }

   private void buildPanel()
   {
      titleLabel.setText("Title:");
      table.setShowGrid(true);
      table.setRequestFocusEnabled(false);
      table.getSelectionModel().addListSelectionListener(this);
      gbc.fill = GridBagConstraints.HORIZONTAL;
      gbc.weighty = 0.0;
      add(titleLabel);
      gbc.fill = GridBagConstraints.BOTH;
      gbc.weighty = 1.0;
      add(scrollPane);
/*
      this.add(titleLabel, new GridBagConstraints(0, 0, 1, 1, 1.0, 0.0
                                                  , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL,
                                                  new Insets(0, 0, 0, 0), 0, 0));
      this.add(scrollPane, new GridBagConstraints(0, 1, 1, 1, 1.0, 1.0
                                                  , GridBagConstraints.CENTER, GridBagConstraints.BOTH,
                                                  new Insets(0, 0, 0, 0), 0, 0));
*/
   }

   public TableCellRenderer getDefaultRenderer(Class columnClass)
   {
	   return table.getDefaultRenderer(columnClass);
   }

   public void setDefaultRenderer(Class columnClass, TableCellRenderer renderer)
   {
	   table.setDefaultRenderer(columnClass, renderer);
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

   public void setTableModel(TableModel tableModel)
   {
	   final TableModel _model = tableModel;
       final StsTablePanel tablePanel = this;
       Runnable runnable = new Runnable()
       {
            public void run()
            {
                table.setModel(_model);
                selectedRows = new byte[_model.getRowCount()];
                for (int i = 0; i < _model.getColumnCount(); i++)
                    table.setDefaultRenderer(table.getColumnClass(i), new StsTableCellRenderer(tablePanel));
            }
       };
       StsToolkit.runWaitOnEventThread(runnable);
   }

   public void addColumns(Object[] headers)
   {
	   int nCols = headers == null ? 0 : headers.length;
	   for (int i = 0; i < nCols; i++)
	   {
		   tableModel.addColumn(headers[i]);
	   }
	   for (int i = 0; i < nCols; i++)
	   {
		   table.setDefaultRenderer(table.getColumnClass(i), new StsTableCellRenderer(this));
	   }
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

   public void changeColumnHeadings(Object[] newHeaders)
   {
      tableModel.setColumnIdentifiers(newHeaders);
   }

   public void removeAllColumns()
   {
      tableModel.setColumnCount(0);
   }

   public void setColumnCount(int nColumns)
   {
      tableModel.setColumnCount(nColumns);
   }

   public void addColumn(Object header)
   {
      tableModel.addColumn(header);
   }

   public void addAccentedColumn(Object header)
   {
      tableModel.addColumn(header);
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
      selectedRows[row] = value;
//        System.out.println("Setting selectedRows[" + row + "]");
   }

   public void setSelectionMode(int selectionMode)
   {
	   table.setSelectionMode(selectionMode);
   }

   public byte getRowType(int row)
   {
	try
	{
		return selectedRows[row];
	}
	catch (Exception ex)
	{
		ex.printStackTrace();
		return (byte)0;
	}
   }

   public void addRows(Object[][] rows)
   {
      int nRows = rows == null ? 0 : rows.length;
      for (int i = 0; i < nRows; i++)
      {
         addRow(rows[i]);
      }
   }

   public void addRow(Object[] row)
   {
      selectedRows = addElementToArray(selectedRows, this.NOT_HIGHLIGHTED);
      tableModel.addRow(row);
   }

   private byte[] addElementToArray(byte[] arrayIn, byte value)
   {
      byte[] arrayOut = null;
      if (arrayIn == null)
         arrayOut = new byte[1];
      else
      {
         arrayOut = new byte[arrayIn.length + 1];
         System.arraycopy(arrayIn, 0, arrayOut, 0, arrayIn.length);
      }
      arrayOut[arrayOut.length - 1] = value;
      return arrayOut;
   }

   public void addData(Object[][] data, Object[] colNames)
   {
      for (int j = 0; j < data.length; j++)
         selectedRows = addElementToArray(selectedRows, this.NOT_HIGHLIGHTED);
      tableModel.setDataVector(data, colNames);
   }

   public int[] getSelectedIndices()
   {
      return table.getSelectedRows();
   }

   public void setColumnSelectable()
   {
      table.setColumnSelectionAllowed(true);
   }

   public void setSelectedIndices(int[] indices)
   {
      DefaultListSelectionModel model = (DefaultListSelectionModel)table.getSelectionModel();
      if (indices == null)
         model.clearSelection();
      else
      {
         for (int i = 0; i < indices.length; i++)
            model.addSelectionInterval(indices[i], indices[i]);
      }
   }

   public void removeAllRows()
   {
      tableModel.setNumRows(0);
      selectedRows = null;
   }

   public synchronized void removeListSelectionListener(ListSelectionListener l)
   {
      if (listSelectionListeners != null && listSelectionListeners.contains(l))
      {
         Vector v = (Vector)listSelectionListeners.clone();
         v.removeElement(l);
         listSelectionListeners = v;
      }
   }

   public synchronized void addListSelectionListener(ListSelectionListener l)
   {
      Vector v = listSelectionListeners == null ? new Vector(2) : (Vector)listSelectionListeners.clone();
      if (!v.contains(l))
      {
         v.addElement(l);
         listSelectionListeners = v;
      }
   }

   protected void fireValueChanged(ListSelectionEvent e)
   {
      if (listSelectionListeners != null)
      {
         Vector listeners = listSelectionListeners;
         int count = listeners.size();
         for (int i = 0; i < count; i++)
            ((ListSelectionListener)listeners.elementAt(i)).valueChanged(e);
      }
   }

   public void valueChanged(ListSelectionEvent e)
   {
      fireValueChanged(e);
   }

   public boolean isSelectedRow(int row)
   {
      for (int i = 0; i < table.getSelectedRowCount(); i++)
      {
         if (row == table.getSelectedRows()[i])
            return true;
      }
      return false;
   }


    public void repaint()
    {
        if(StsPostStackTraceDefinitionPanel.debug) StsException.systemDebug(this, "repaint", " time:" + System.currentTimeMillis());
        super.repaint();
    }

   static public void main(String[] args)
   {
   /*
       Object[] rowObjects;
       int nObjects = 3;
       rowObjects = new Object[nObjects];
       String[] columnNames = new String[] { "xMin", "xMax", "yMin", "yMax" };
       for(int n = 0; n < nObjects; n++)
       {
           rowObjects[n] = new StsBoundingBox(n+0.0f, n+0.1f, n+0.2f, n+0.3f);
       }
       StsObjectTableModel2 tableModel = new StsObjectTableModel2(rowObjects, columnNames);
       StsTablePanel tablePanel = new StsTablePanel(tableModel);
       StsToolkit.createDialog(tablePanel);
   */
      JDialog d = new JDialog((Frame)null, "Table Panel Test", true);
      StsTablePanel panel = new StsTablePanel();

      d.getContentPane().add(panel);

      panel.setTitle("title of panel .....");
      Object[] colnames =
         {"One", "Two", "Three"};
      for (int i = 0; i < 3; i++)
         panel.addColumn(colnames[i]);

      Object[] row = new Object[3];
      for (int j = 0; j < 10; j++)
      {
         for (int i = 0; i < 3; i++)
            row[i] = new String("(" + i + "," + j + ")");
         panel.addRow(row);
      }

      d.pack();
      d.setVisible(true);
   }

   class StsTableCellRenderer extends DefaultTableCellRenderer
   {
      StsTablePanel panel = null;

      Color bg_highlighted = new Color(255, 255, 165);
      Color fg_highlighted = Color.BLACK;
      Color bg_not_highlighted = Color.WHITE;
      Color fg_not_highlighted = Color.BLACK;
      Color bg_not_editable = new Color(220, 220, 220);
      Color fg_not_editable = Color.BLACK;
      Color bg_selected = Color.BLUE;
      Color fg_selected = Color.WHITE;

      public StsTableCellRenderer(StsTablePanel panel)
      {
         this.panel = panel;
      }

      // method to override - returns cell renderer component
      public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,
         int row, int column)
      {
         // let the default renderer prepare the component for us
         Component comp = super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);

         // now get the current font used for this cell
         Font font = comp.getFont();
         if (panel.isSelectedRow(row))
         {
            comp.setBackground(bg_selected);
            comp.setForeground(fg_selected);
         }
         else if (panel.getRowType(row) == StsTablePanelNew.HIGHLIGHTED)
         {
            comp.setBackground(bg_highlighted);
            comp.setForeground(fg_highlighted);
         }
         else if (panel.getRowType(row) == StsTablePanelNew.NOT_EDITABLE)
         {
            comp.setBackground(bg_not_editable);
            comp.setForeground(fg_not_editable);
         }
         else if (panel.getRowType(row) == StsTablePanelNew.NOT_HIGHLIGHTED)
         {
            comp.setBackground(bg_not_highlighted);
            comp.setForeground(fg_not_highlighted);
         }
         return comp;
      }
   }
}
