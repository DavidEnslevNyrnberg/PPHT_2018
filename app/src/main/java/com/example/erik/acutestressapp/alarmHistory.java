package com.example.erik.acutestressapp;

import android.os.Bundle;
import android.support.v7.app.AppCompatActivity;
import android.widget.ArrayAdapter;
import android.widget.ListView;

public class alarmHistory extends AppCompatActivity {

	// Creates the layout
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.alarm_history);
		final ListView listview = findViewById(R.id.listview);

		// At this time the Alarm history is just a place holder place where there's a list of
		// strings with time stamps. In future when a SVM alert is set it will append the time
		// stamp in to a array of strings, and be displayed here.

		String[] values = new String[] { "08.04.18 - 12:25", "08.04.18 - 14:35", "08.04.18 - 16:10",
				"10.04.18 - 11:05", "01.05.18 - 22:10"};
				// List the values in the String array to a
		ArrayAdapter<String> adapter = new ArrayAdapter<String>(this,R.layout.list_white_text,
				R.id.list_content, values);
		listview.setAdapter(adapter);


	}
}