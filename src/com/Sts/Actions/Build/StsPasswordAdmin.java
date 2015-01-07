
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.Utilities.*;

import javax.crypto.*;
import javax.swing.*;
import java.awt.*;
import java.io.*;

public class StsPasswordAdmin extends JPanel
{	    
	JTextArea textArea = new JTextArea("Generating key file and password...");
 	public StsPasswordAdmin() 
 	{
        super(new BorderLayout());
 		this.add(textArea);
 		
 		SecretKey key = null;
		String homeDirectory = System.getProperty("user.home");		
 		try 
 		{
 			key = createSecretKey();
 		}
 		catch (Exception e) 
 		{
 			textArea.setText("Problem creating and saving password key" + e.toString());
 			e.printStackTrace();
 			return;
 		}
 		
		// Output key to file - s2ssystems.key
 		FileOutputStream fos = null;
 	    ObjectOutputStream out = null;
 	    try
 	    {
 	    	fos = new FileOutputStream(homeDirectory + "/s2ssystems.key");
 	    	out = new ObjectOutputStream(fos);
 	    	out.writeObject(key);
 	    	out.writeLong(System.currentTimeMillis());
 	    	out.close();
 	    }
 	    catch(IOException ex)
 	    {
 			textArea.setText("Problem saving password key to file:\n\n" + ex.toString());
 			ex.printStackTrace();
 			return;
 	    }
		
 		// Produce password (key and password must both be provided to the user)
 	    StsPassword passwordTool = null;
 	    try
 	    {
 	    	passwordTool = new StsPassword(key);
 	    }
 	    catch(Exception ex)
 	    {
 	    	textArea.setText("Problem creating password object" + ex.toString());
 	    	ex.printStackTrace();
 	    	return;
 	    } 
 		long userDate = System.currentTimeMillis();
 	    String encryptUser = passwordTool.encrypt("S2S" + Long.toString(userDate)); 			
 		
 	    // Inform user of the new password and location of the key file.
 	    textArea.setText("Key file is named:\n\n     " +
 	    		"s2ssystems.key\n\n" +
 	    		"It is located at:\n\n     " + homeDirectory +
 	    		"\n\nThe encrypted password is:\n\n     " + encryptUser);
 	}
 	
    private static void createAndShowGUI() 
    {
        //Create and set up the window.
        JFrame frame = new JFrame("Key and Password Generation Tool");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        //Create and set up the content pane.
        JComponent newContentPane = new StsPasswordAdmin();
        newContentPane.setOpaque(true); //content panes must be opaque
        frame.setContentPane(newContentPane);

        //Display the window.
        frame.pack();
        frame.setVisible(true);
    }
    
 	private static SecretKey createSecretKey()
 	{
 		try
 		{
 			KeyGenerator keyGen = KeyGenerator.getInstance("DES");
 			return keyGen.generateKey();
 		}
 		catch(Exception ex)
 		{
 			System.out.println("Unable to get instance of key generator....");
 			return null;
 		}
 	}
 	
    public static void main(String[] args)
    {
        StsToolkit.runLaterOnEventThread(new Runnable() { public void run() { createAndShowGUI(); }});
    } 	
} 

