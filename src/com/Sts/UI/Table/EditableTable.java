package com.Sts.UI.Table;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 14, 2007
 * Time: 3:39:23 PM
 * To change this template use File | Settings | File Templates.
 */

import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.table.*;
import java.lang.reflect.*;
import java.util.*;

public class EditableTable
{

    public static void main(String[] a)
    {
        /*
        Runnable runnable = new Runnable()
        {
            public void run()
            {
        */
                JFrame frame = new JFrame();
                frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

                DataEntry[] dataEntries =  new DataEntry[]
                {
                    new DataEntry("Saravan", "Pantham", new Integer(50), "B", new Boolean(false)),
                    new DataEntry("Eric", "", new Integer(180), "O", new Boolean(true)),
                    new DataEntry("John", "", new Integer(120), "AB", new Boolean(false)),
                    new DataEntry("Mathew", "", new Integer(140), "A", new Boolean(true))
                };

                ArrayList entryList = new ArrayList();
                for(int n = 0; n < dataEntries.length; n++)
                    entryList.add(dataEntries[n]);

                EditableTableModel model = new EditableTableModel(entryList);
                JTable table = new JTable(model);
                model.setTable(table);
                table.createDefaultColumnsFromModel();

                String[] bloodGroups = {"A", "B", "AB", "O"};
                JComboBox comboBox = new JComboBox(bloodGroups);
                table.getColumnModel().getColumn(3).setCellEditor(new DefaultCellEditor(comboBox));

                frame.add(new JScrollPane(table));
                frame.setSize(300, 200);
                frame.setVisible(true);
                table.repaint();
                StsToolkit.sleep(2000);

                DataEntry newEntry =  new DataEntry("Tom", "", new Integer(220), "A", new Boolean(true));
                model.addRow(newEntry);
                table.repaint();
        /*
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
        */
    }
}

class DataEntry implements StsTableRowInterface
{
    public String firstName;
    public String lastName;
    public int weight;
    public String bloodGroup;
    public boolean old;

    static String[] columnTitles = {"First Name", "Last Name", "Weight (lb)", "Blood Group", "Age>20yrs"};
    static String[] columnNames = {"firstName", "lastName", "weight", "bloodGroup", "old"};

    public DataEntry(String firstName, String lastName, int weight, String bloodGroup, boolean old)
    {
        this.firstName = firstName;
        this.lastName = lastName;
        this.weight = weight;
        this.bloodGroup = bloodGroup;
        this.old = old;
    }

    public String[] getColumnTitles() { return columnTitles; }
    public String[] getColumnNames() { return columnNames; }
}

interface StsTableRowInterface
{
    public String[] getColumnNames();
    public String[] getColumnTitles();
}

class EditableTableModel extends AbstractTableModel
{
    JTable table;
    String[] columnTitles;
    String[] columnNames;
    ArrayList rowObjects;
    Field[] columnFields;
    Class rowClass;
    int nRows;
    int nCols;
    int nDeleteButtonColumn;

    public EditableTableModel(ArrayList rowObjects)
    {
        StsTableRowInterface rowObject = (StsTableRowInterface)rowObjects.get(0);
        this.rowObjects = rowObjects;
        columnTitles = rowObject.getColumnTitles();
        columnNames = rowObject.getColumnNames();
        this.rowClass = rowObject.getClass();
        nRows = rowObjects.size();
        nCols = columnTitles.length;
        initializeColumnFields();
        addDeleteButtons();
    }

    public void setTable(JTable table)
    {
        this.table = table;
    }

    private void initializeColumnFields()
    {
        columnFields = new Field[nCols];
        for(int n = 0; n < nCols; n++)
        {
            try
            {
                columnFields[n] = rowClass.getDeclaredField(columnNames[n]);
            }
            catch(Exception e)
            {
                StsException.systemError(this, "constructor", "Failed to find field " + columnNames[n] + " in class " + rowClass.getName());
            }
        }
    }

    public void addDeleteButtons()
    {
        nDeleteButtonColumn = nCols;
        nCols++;
        columnTitles = (String[])StsMath.arrayAddElement(columnNames, "delete");
    }

    public int getRowCount()
    {
        return rowObjects.size();
    }

    public int getColumnCount()
    {
        return columnTitles.length;
    }

    public Object getValueAt(int row, int col)
    {
        try
        {
            if(col == nDeleteButtonColumn)
                return new Boolean(false);
            else
                return columnFields[col].get(rowObjects.get(row));
        }
        catch(Exception e)
        {
            StsException.systemError(this, "getValueAt", "Failed to find value at row " + row + " col " + col);
            return null;
        }
    }

    public String getColumnName(int col)
    {
        return columnTitles[col];
    }

    public Class getColumnClass(int col)
    {
        if(col == nDeleteButtonColumn)
            return Boolean.class;
        else
            return columnFields[col].getDeclaringClass();
//        return getValueAt(0, col).getClass();
    }

    public boolean isCellEditable(int row, int col)
    {
        return true;
    }

    public void setValueAt(Object value, int row, int col)
    {
       try
        {
           if(col == nDeleteButtonColumn)
                deleteRow(row);
            else
               columnFields[col].set(rowObjects.get(row), value);
        }
        catch(Exception e)
        {
        }
    }

    public void addRow(Object rowObject)
    {
        if(!rowClass.isInstance(rowObject))
        {
            System.out.println("Can't add rowObject of class " + StsToolkit.getSimpleClassname(rowObject) + ". All objects must be of class " + rowClass.getName());
        }
        rowObjects.add(rowObject);
        nRows++;
        table.repaint();
    }

    private void deleteRow(int row)
    {
        rowObjects.remove(row);
        nRows--;
        System.out.println("Delete row: " + row);
        table.repaint();
    }
}
