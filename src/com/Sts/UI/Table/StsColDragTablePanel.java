package com.Sts.UI.Table;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.table.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jan 21, 2010
 * Time: 3:14:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class StsColDragTablePanel extends StsJPanel
{
    int nCols;
    DefaultTableModel tableModel;
    JTable table;
    static String dragSymbol = new String("<--->");

   public StsColDragTablePanel(String[] colNames)
   {
        nCols = colNames.length;
        GridLayout gridLayout = new GridLayout(1, nCols);
        JPanel headerPanel = new JPanel(gridLayout);

        Object[] dragSymbols = new String[nCols];
        for(int n = 0; n < nCols; n++)
        {
            JLabel colHeader = new JLabel(colNames[n], SwingConstants.CENTER);
            headerPanel.add(colHeader);
            dragSymbols[n] = dragSymbol;
        }
        tableModel = new DefaultTableModel(dragSymbols, 0);
        table = new JTable(tableModel);
        table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
        gbc.fill = GridBagConstraints.BOTH;
        gbc.weighty = 0.0;
        add(headerPanel);

        JScrollPane scrollPane = new JScrollPane(table);
        gbc.weighty = 1.0;
        add(scrollPane);
   }

   public void addRow(Object[] row)
   {
       tableModel.addRow(row);
   }

   static public void main(String[] args)
   {
      Runnable runnable = new Runnable()
      {
          public void run()
          {
              int ncols = 3;
              int nrows = 10;
              String[] colNames = new String[] { "one", "Two", "Three"};
              StsColDragTablePanel tablePanel = new StsColDragTablePanel(colNames);

              Object[] row = new Object[ncols];
              for (int i = 0; i < nrows; i++)
              {
                 for (int j = 0; j < ncols; j++)
                    row[j] = new String("(" + i + "," + j + ")");
                 tablePanel.addRow(row);
              }
              StsToolkit.createDialog(tablePanel);
          }
      };
      StsToolkit.runWaitOnEventThread(runnable);
   }
}
