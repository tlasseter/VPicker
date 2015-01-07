
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Build;

import com.Sts.UI.*;

import javax.crypto.*;
import java.awt.*;
import java.io.*;

public class StsPassword
{
 	Cipher ecipher;
 	Cipher dcipher;

 	public StsPassword(SecretKey key) {
 		try
 		{
 			ecipher = Cipher.getInstance("DES");
 			dcipher = Cipher.getInstance("DES");
 			ecipher.init(Cipher.ENCRYPT_MODE, key);
 			dcipher.init(Cipher.DECRYPT_MODE, key);

 		} catch (javax.crypto.NoSuchPaddingException e) {
 		} catch (java.security.NoSuchAlgorithmException e) {
 		} catch (java.security.InvalidKeyException e) {
 		}
 	}

 	public String encrypt(String str) {
 		try {
 			// Encode the string into bytes using utf-8
 			byte[] utf8 = str.getBytes("UTF8");

 			// Encrypt
 			byte[] enc = ecipher.doFinal(utf8);

 			// Encode bytes to base64 to get a string
 			return new sun.misc.BASE64Encoder().encode(enc);
 		} catch (javax.crypto.BadPaddingException e) {
 		} catch (IllegalBlockSizeException e) {
 		} catch (UnsupportedEncodingException e) {
 		} catch (java.io.IOException e) {
 		}
 		return null;
 	}

 	public String decrypt(String str) {
 		try {
 			// Decode base64 to get bytes
 			byte[] dec = new sun.misc.BASE64Decoder().decodeBuffer(str);

 			// Decrypt
 			byte[] utf8 = dcipher.doFinal(dec);

 			// Decode using utf-8
 			return new String(utf8, "UTF8");
 		} catch (javax.crypto.BadPaddingException e) {
 		} catch (IllegalBlockSizeException e) {
 		} catch (UnsupportedEncodingException e) {
 		} catch (java.io.IOException e) {
 		}
 		return null;
 	}

 	public boolean writeLastUsed(SecretKey key)
 	{
		// Read the key from the users home directory
	    String filename = "s2ssystems.key";
		String homeDirectory = System.getProperty("user.home");

 		FileOutputStream fos = null;
 	    ObjectOutputStream out = null;
 	    try
 	    {
 	    	fos = new FileOutputStream(homeDirectory + "/" + filename);
 	    	out = new ObjectOutputStream(fos);
 	    	out.writeObject(key);
 	    	out.writeLong(System.currentTimeMillis());
 	    	out.close();
 	    }
	    catch(Exception ex)
        {
	    	new StsMessage(new Frame(), StsMessage.WARNING, "Unable to write last used date to key file");
 		    ex.printStackTrace();
 		    return false;
 	    }
 		return true;
 	}

 	public long readLastUsed()
 	{
 		long lastUsedOn = 0L;

		// Read the key from the users home directory
	    String filename = "s2ssystems.key";
		String homeDirectory = System.getProperty("user.home");

	    FileInputStream fis = null;
	    ObjectInputStream in = null;
	    try
	    {
	    	fis = new FileInputStream(homeDirectory + "/" + filename);
	    	in = new ObjectInputStream(fis);
	    	in.readObject();
	    	lastUsedOn = in.readLong();
	    }
	    catch(Exception ex)
        {
	    	new StsMessage(new Frame(), StsMessage.WARNING, "Unable to read last used date from key file");
 		    ex.printStackTrace();
 		    return 0L;
 	    }
	    return lastUsedOn;
 	}

 	public static SecretKey readSecretKey()
 	{
 		SecretKey key = null;

		// Read the key from the users home directory
	    String filename = "s2ssystems.key";
		String homeDirectory = System.getProperty("user.home");

	    FileInputStream fis = null;
	    ObjectInputStream in = null;
	    try
	    {
	    	fis = new FileInputStream(homeDirectory + "/" + filename);
	    	in = new ObjectInputStream(fis);
	    	key = (SecretKey)in.readObject();
	    	in.close();
	    }
	    catch(Exception ex)
        {
	    	new StsMessage(new Frame(), StsMessage.WARNING, "Unable to read key file from home directory");
 		    ex.printStackTrace();
 		    key = null;
 	    }

 		// Key file is not found
 		if(key == null)
 		{
 			if(!StsYesNoDialog.questionValue(new Frame(), "Unable to locate or read key file:\n\n    " +
 					homeDirectory + "/" + filename +
 					"\n\nWould you like to select a different file?"))
 				return null;
 			else
 			// Show file selection dialog
 				System.out.println("Add file selection dialog...");
 		}

 		return key;
 	}
}

