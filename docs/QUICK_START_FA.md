# راهنمای کوتاه اجرا

۱. فایل ZIP را استخراج کن و وارد پوشهٔ `PanCancer-miRNA-CYT` شو.

۲. برای بررسی بدون اجرای تحلیل و بدون دریافت داده‌های حجیم:

```bash
Rscript scripts/run.R --dry-run
Rscript scripts/install_dependencies.R --tests-only
Rscript tests/run_tests.R
```

۳. اگر بسته‌ها موجود نیستند:

```bash
Rscript scripts/install_dependencies.R --install
```

۴. داده‌ها را از پوشهٔ پروژهٔ قبلی منتقل کن؛ مسیر نمونه را عوض کن:

```bash
Rscript scripts/prepare_inputs.R --from=/path/to/old/MicroRNA
Rscript scripts/check_environment.R
```

۵. اجرای کامل:

```bash
Rscript scripts/run.R --stages=1-27
```

اگر مثلاً مرحلهٔ ۲۱ متوقف شد، پس از رفع علت، همان مرحله و ادامه را اجرا کن:

```bash
Rscript scripts/run.R --stages=21-27
```

مرحلهٔ قدیمی ۱۳ حذف شده است؛ مرحلهٔ جدید ۱۳ متاآنالیز است. خوشه‌بندی سنگین
اکنون مرحلهٔ ۱۷ است. انتقال نتیجهٔ قبلی آن به فایل‌های کامل همان تحلیل نیاز
دارد؛ فقط گزارش یا تصویر برای این کار کافی نیست. دستور انتقال در README است.

برای RStudio فایل پروژه را باز کن و بنویس:

```r
source("scripts/run.R")
run_pipeline(stages = 1:27)
```

خروجی‌ها در `results/` و شکل‌ها و جدول‌های نهایی در `results/manuscript/`
قرار می‌گیرند. برای GitHub محتوای همین پوشهٔ کد را استفاده کن؛ داده‌های حجیم و
خروجی‌های محلی توسط `.gitignore` کنار گذاشته شده‌اند. اجرای کامل دوباره روی
داده‌های واقعی در محیط آماده‌سازی این بسته انجام نشده است؛ جزئیات آزمون‌ها
در `VALIDATION.md` ثبت شده‌اند.
