package com.Sts.UI;

import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsNameSelectionDialog extends JDialog
{
   public static boolean ALLOW_CANCEL = true;
   public static boolean QUESTION_OVERWRITE = true;

   private StsJPanel panel = StsJPanel.addInsets();
   private JPanel buttonPanel = new JPanel();
   private JLabel help = new JLabel("Select a name or enter a new one.");
   private JLabel label = new JLabel();
   private JComboBox comboBox = new JComboBox();
   private JButton cancelButton = null;
   private JButton okayButton = new JButton("OK");

   private boolean questionOverwrite;
   private boolean wasCanceled = false;

   public StsNameSelectionDialog(Frame owner, String title, String listName, String[] names, String selectedName,
                                 boolean questionOverwrite, boolean allowCancel)
   {
      super(owner, title, true);
      this.questionOverwrite = questionOverwrite;

      try
      {
         setDefaultCloseOperation(DISPOSE_ON_CLOSE);
         jbInit(allowCancel);

         label.setText(listName);
         comboBox.setEditable(true);
         DefaultComboBoxModel comboModel = new DefaultComboBoxModel(names);
         comboBox.setModel(comboModel);
         comboModel.setSelectedItem(selectedName);
		 StsToolkit.centerComponentOnScreen(this);

         pack();
      }
      catch (Exception exception)
      {
         exception.printStackTrace();
      }
   }

   public String getSelectedName()
   {
      return (String)comboBox.getSelectedItem();
   }

   public boolean wasCanceled()
   {
      return wasCanceled;
   }

   private void jbInit(boolean allowCancel)
      throws Exception
   {
      panel.addEndRow(help, 4, 0.0);
      panel.addToRow(label, 1, 0.0);
      panel.addEndRow(comboBox, 3, 1.0);
      buttonPanel.setLayout(new FlowLayout());
      panel.addEndRow(buttonPanel, 4, 1.0);
      buttonPanel.add(okayButton);
      if (allowCancel)
      {
         cancelButton = new JButton("Cancel");
         buttonPanel.add(cancelButton);
         cancelButton.addActionListener(new java.awt.event.ActionListener()
         {
            public void actionPerformed(ActionEvent e)
            {
               cancelButtonAction();
            }
         });
      }
      okayButton.addActionListener(new java.awt.event.ActionListener()
      {
         public void actionPerformed(ActionEvent e)
         {
            okayButtonAction();
         }
      });

      getContentPane().add(panel);
   }

   private void okayButtonAction()
   {
      if (questionOverwrite && comboBox.getSelectedIndex() > -1)
      {
         boolean overwrite = StsYesNoDialog.questionValue(this.getOwner(),
            "Name exists: do you wish to overwrite it?");
         if (!overwrite)
            return;
      }
      setVisible(false);
   }

   private void cancelButtonAction()
   {
      wasCanceled = true;
      setVisible(false);
   }

   public static void main(String[] args)
   {
      String[] items = new String[]
         {"One", "Two"};
      StsNameSelectionDialog d = new StsNameSelectionDialog(null, "StsComboBoxDialog test", "list of items:", items,
         "Two", true, true);
      d.setVisible(true);
      System.exit(0);
   }
}
