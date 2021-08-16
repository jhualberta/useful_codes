from bs4 import BeautifulSoup
import requests

# 设置http请求头伪装成浏览器
headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/61.0.3163.100 Safari/537.36",
    "Connection": "keep-alive",
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8",
    "Accept-Language": "zh-CN,zh;q=0.8"}

# requests获取博客页面html文本
url = "https://www.douban.com/people/sciama/notes"
#https://www.douban.com/people/sciama/notes
#https://www.douban.com/people/sciama/statuses
data = {
    'username':'',   
    'password':'',
}
#r = requests.get(url, data, headers=send_headers)
session = requests.Session()
session.post(url, headers = headers,data = data)
r = session.get(url, headers = headers, data = data)
r.encoding = "utf-8"
html = r.text

# 将获取到的html送入bs4进行解析
soup = BeautifulSoup(html, "html.parser")   # 获得解析后对象
content = soup.find("div", id="content")    # 找到id是content的div
print(content)
# 找到这个div中所有 class 是 article-item-box csdn-tracking-statistics 的div
artlist = content.find_all("div", attrs={"class":"note-header-container"})
# 遍历每个div 输出内容 以html形式输出
for div in artlist:
    a = div.h3.a
    print("<a href='" + a["href"] + "'> " + a["title"])
