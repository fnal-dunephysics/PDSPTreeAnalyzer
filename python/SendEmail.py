import smtplib
from email.mime.text import MIMEText

def SendEmail(From,To,Subject,Content):

  msg = MIMEText(Content)

  msg['Subject'] = Subject
  msg['From'] = From
  msg['To'] = To

  s = smtplib.SMTP('localhost')
  s.sendmail(From, [To], msg.as_string())
  s.quit()
