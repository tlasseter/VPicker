package com.Sts.UI.Table;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import javax.swing.event.*;
import javax.swing.table.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;

public class EditableBeanTable
{

    public static void main(String[] a)
    {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        BeanDataEntry[] dataEntries =  new BeanDataEntry[]
        {
            new BeanDataEntry("Saravan", "Pantham", new Integer(50), "B", new Boolean(false)),
            new BeanDataEntry("Eric", "", new Integer(180), "O", new Boolean(true)),
            new BeanDataEntry("John", "", new Integer(120), "AB", new Boolean(false)),
            new BeanDataEntry("Mathew", "", new Integer(140), "A", new Boolean(true))
        };

        EditableBeanTableModel model = new EditableBeanTableModel(dataEntries);
        BeanTable table = new BeanTable(model);
        frame.add(new JScrollPane(table));
        frame.setSize(300, 200);
        frame.setVisible(true);
    }
}

class BeanTable extends JTable
{
    BeanTable(EditableBeanTableModel model)
    {
        super(model);
    }

    public void createDefaultColumnsFromModel()
    {
        EditableBeanTableModel m = (EditableBeanTableModel)getModel();
        if (m != null) 
        {
            // Remove any current columns
            TableColumnModel cm = getColumnModel();
            while (cm.getColumnCount() > 0)
            {
                cm.removeColumn(cm.getColumn(0));
	        }

            // Create new columns from the data model info
            StsFieldBean[] columnBeans = m.columnBeans;
            for (int i = 0; i < m.getColumnCount(); i++)
            {
                BeanCellRendererEditor beanCell = new BeanCellRendererEditor(columnBeans[i]);
                addColumn(new TableColumn(i, 10, beanCell, beanCell));
//                BeanCellEditor cellEditor = BeanCellEditor.constructor(columnBeans[i]);
//                addColumn(new TableColumn(i, 10, columnBeans[i], cellEditor));
            }
        }
    }
}

class BeanCellEditor extends DefaultCellEditor
{
    BeanCellEditor(StsTextFieldBean stringBean)
    {
        super(stringBean.getTextField());
    }

    BeanCellEditor(StsComboBoxFieldBean comboBoxBean)
    {
        super(comboBoxBean.getComboBox());
        comboBoxBean.getComboBox().putClientProperty("JComboBox.isTableCellEditor", Boolean.TRUE);
    }

    BeanCellEditor(StsBooleanFieldBean booleanBean)
    {
        super(booleanBean.getCheckBox());
    }

    static BeanCellEditor constructor(StsFieldBean fieldBean)
    {
        if(fieldBean instanceof StsTextFieldBean)
            return new BeanCellEditor((StsTextFieldBean)fieldBean);
        else if(fieldBean instanceof StsComboBoxFieldBean)
            return new BeanCellEditor((StsComboBoxFieldBean)fieldBean);
        else if(fieldBean instanceof StsBooleanFieldBean)
            return new BeanCellEditor((StsBooleanFieldBean)fieldBean);
        else
        {
            StsException.systemError("BeanCellEditor constructor failed. No constructor for bean of type " + StsToolkit.getSimpleClassname(fieldBean));
            return null;
        }
    }
}

class BeanCellRendererEditor extends AbstractCellEditor implements TableCellEditor, TableCellRenderer
{
    StsFieldBean bean;
    int clickCountToStart = 0;

    BeanCellRendererEditor(StsFieldBean bean)
    {
        this.bean = bean;
    }

    public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus,  int row, int column)
    {
        Object valueObject = table.getModel().getValueAt(row, column);
        bean.setValueObject(valueObject);
        return bean.getMainComponent();
    }

   public Component getTableCellEditorComponent(JTable table, Object value, boolean isSelected, int row, int column)
   {
       Object valueObject = table.getModel().getValueAt(row, column);
       bean.setValueObject(valueObject);
       return bean.getMainComponent();
   }

    public Object getCellEditorValue()
    {
        return bean.getValueObject();
    }

    public boolean isCellEditable(EventObject anEvent)
    {
        System.out.println("isCellEditable event: " + anEvent.toString());
        if (anEvent instanceof MouseEvent)
        {
            return ((MouseEvent)anEvent).getClickCount() >= clickCountToStart;
        }
        return bean.getEditable();
    }

    public boolean shouldSelectCell(EventObject anEvent)
    {
        System.out.println("shouldSelectCell event: " + anEvent.toString());
        if (anEvent instanceof MouseEvent)
        {
            MouseEvent e = (MouseEvent)anEvent;
            return e.getID() != MouseEvent.MOUSE_DRAGGED;
        }
        return bean.getEditable();
    }

    public boolean stopCellEditing()
    {
        System.out.println("stopCellEditing");
//        fireEditingStopped();
	    return true;
    }

    public void  cancelCellEditing()
    {
        System.out.println("cancelCellEditing");
//        fireEditingCanceled();
    }
}

class BeanDataEntry implements StsBeanTableRowInterface
{
    private String firstName;
    private String lastName;
    private int weight;
    private String bloodGroup;
    private boolean old;

    static String[] bloodGroups = new String[] { "A", "B", "O", "AB" };

    static StsFieldBean[] entryBeans = new StsFieldBean[]
    {
            new StsStringFieldBean(BeanDataEntry.class, "firstName", true, null),
            new StsStringFieldBean(BeanDataEntry.class, "lastName", true, null),
            new StsIntFieldBean(BeanDataEntry.class, "weight", true, null),
            new StsComboBoxFieldBean(BeanDataEntry.class, "bloodGroup", null, bloodGroups),
            new StsBooleanFieldBean(BeanDataEntry.class, "old", true, null)
    };
    static String[] columnTitles = {"First Name", "Last Name", "Weight (lb)", "Blood Group", "Age>20yrs"};

    public BeanDataEntry(String firstName, String lastName, int weight, String bloodGroup, boolean old)
    {
        this.firstName = firstName;
        this.lastName = lastName;
        this.weight = weight;
        this.bloodGroup = bloodGroup;
        this.old = old;
    }

    public String[] getColumnTitles() { return columnTitles; }
    public StsFieldBean[] getColumnBeans() { return entryBeans; }

    public String getFirstName()
    {
        return firstName;
    }

    public void setFirstName(String firstName)
    {
        this.firstName = firstName;
    }

    public String getLastName()
    {
        return lastName;
    }

    public void setLastName(String lastName)
    {
        this.lastName = lastName;
    }

    public int getWeight()
    {
        return weight;
    }

    public void setWeight(int weight)
    {
        this.weight = weight;
    }

    public String getBloodGroup()
    {
        return bloodGroup;
    }

    public void setBloodGroup(String bloodGroup)
    {
        this.bloodGroup = bloodGroup;
    }

    public boolean getOld()
    {
        return old;
    }

    public void setOld(boolean old)
    {
        this.old = old;
    }
}

interface StsBeanTableRowInterface
{
    public StsFieldBean[] getColumnBeans();
    public String[] getColumnTitles();
}

class EditableBeanTableModel implements TableModel
{
    String[] columnTitles;
    StsFieldBean[] columnBeans;
    Object[] rowObjects;
    Class rowClass;
    int nRows;
    int nCols;
    int nDeleteButtonColumn = -1;
    boolean[] delete;

    public EditableBeanTableModel(StsBeanTableRowInterface[] rowObjects)
    {
        StsBeanTableRowInterface rowObject = rowObjects[0];
        this.rowObjects = rowObjects;
        columnTitles = rowObject.getColumnTitles();
        columnBeans = rowObject.getColumnBeans();
        this.rowClass = rowObject.getClass();
        nRows = rowObjects.length;
        nCols = columnTitles.length;
//        addDeleteButtons();
    }
/*
    private void initializeColumnFields()
    {
        columnFields = new Field[nCols];
        for(int n = 0; n < nCols; n++)
        {
            try
            {
                columnFields[n] = rowClass.getField(columnBeans[n]);
            }
            catch(Exception e)
            {
                StsException.systemError(this, "constructor", "Failed to find field " + columnBeans[n] + " in class " + rowClass.getName());
            }
        }
    }
*/
    public void addDeleteButtons()
    {
        nDeleteButtonColumn = nCols;
        nCols++;
        columnTitles = (String[])StsMath.arrayAddElement(columnBeans, "delete");
        delete = new boolean[nRows];
    }

    public int getRowCount()
    {
        return rowObjects.length;
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
                return delete[row];
            else
                return columnBeans[col].getValueFromPanelObject(rowObjects[row]);
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
        return columnBeans[col].getClass();
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
           {
               delete[row] = ((Boolean)value).booleanValue();
               if(delete[row]) deleteRow(row);
           }
            else
           {
               columnBeans[col].setValueInPanelObject(rowObjects[row], value);
           }
        }
        catch(Exception e)
        {
        }
    }

    public void addTableModelListener(TableModelListener l)
    {

    }

    public void removeTableModelListener(TableModelListener l)
    {

    }

    private void deleteRow(int row)
    {
//        System.out.println("Delete row: " + row);
    }
}
